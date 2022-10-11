module app;

import core.atomic;
import core.memory;

import std.algorithm;
import std.complex;
import std.conv;
import std.exception;
import std.file;
import std.format;
import std.json;
import std.math;
import std.meta;
import std.parallelism;
import std.path;
import std.random;
import std.range;
import std.range;
import std.stdio;
import std.typecons;

import dffdd.utils.unit;
import dffdd.blockdiagram.noise;
import dffdd.dpd.polynomial : DPDMode;

import models;
import simmain;

import msgpack;

import tuthpc.tasklist;
import tuthpc.cluster;
import tuthpc.tasklist;

extern(C) void openblas_set_num_threads(int num_threads);
extern(C) int openblas_get_num_threads();


enum NO_NEED_CSV_FILES = true;
immutable string ALL_RESULT_FILENAME = "all_result.bin";
immutable string DUMPED_RESULT_LIST_FILENAME = "dumped_results.bin";

void writeJSONData(string path, JSONValue jv)
{
    std.file.write(path, msgpack.pack(msgpack.fromJSONValue(jv)));
    std.file.write(path.withExtension(".json"), jv.toString(JSONOptions.specialFloatLiterals));
}


JSONValue readJSONData(string path)
{
    return msgpack.toJSONValue(msgpack.unpack(cast(ubyte[]) std.file.read(path)));
}



shared static this()
{
    import tuthpc.cluster;

    if(auto inst = cast(KyotoBInfo)ClusterInfo.currInstance)
    {
        if(inst.maxPPN == 36)
            openblas_set_num_threads(2);
        else
            openblas_set_num_threads(1);
    }
    else
        openblas_set_num_threads(1);
}



struct ModelSeed
{
    string cancellerType;

    bool linearMode = false;
    bool useIQImbalance = false;
    bool usePhaseNoise = false;
    bool useDPD = true;
    // bool syncPhaseNoise = true;

    /* IQ Mixer */
    Gain txIRR = 25.dB;
    Gain rxIRR = 25.dB;
    double pn3dBBWHz = 0.1;

    /* PA */
    string[2] paModelName;
    Voltage[2] txPower = 23.dBm;       // 送信電力
    Voltage[2] paVsat = 30.dBm;        // 出力飽和電圧
    Gain[2] paGain = 30.dB;
    double[2] paSmoothFactor = 3;

    /* LNA */
    string[2] lnaModelName;
    double[2] lnaSmoothFactor = 3;
    Gain[2] lnaGain = 20.dB;         // LNA利得 
    Voltage[2] lnaVsatIn = (-6).dBm; // 入力飽和電圧

    /* Other */
    Nullable!Gain SNR;
    Nullable!Gain DesiredLOSS;
    Nullable!Gain INR;
    Nullable!Gain TXRXISO;
    uint numOfTrainingSymbols = 20;
    Gain gamma = 0.dB;

    /* DPD Param */
    uint orderOfDPD = 3;
    DPDMode dpdMode = DPDMode.Linearity;

    /* NLMS/RLS Params */
    real nlmsMu;
    real rlsDelta;
    real rlsLambda;

    /* Basis function selection */
    Gain bfsNoiseMargin = (-4).dB;
    uint bfsNumOfEst = 2;

    /* frequency domain iterative */
    uint iterNumOfIteration = 3;
    uint iterNumOfNewton = 1;
    uint iterNumOfSCForEstNL = 8;
    string iterEstOrder = "IHD";
    Flag!"use3rdSidelobe" iterUse3rdSidelobe = Yes.use3rdSidelobe;

    /* measure only desired signal */
    bool onlyDesired = false;

    size_t numofTapsOfSIChannel = 64;
    size_t numOfTapsOfDesiredChannel = 64;

    bool outputBER = false;
}



void mainJob()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    auto setRLSLMSParam(ref ModelSeed model)
    {
        if(model.cancellerType.startsWith("PH")) {
            model.rlsDelta = 0.1;
            model.rlsLambda = 1;
            model.nlmsMu = 0.5;
        } else {
            model.rlsDelta = 4E-7;
            model.rlsLambda = 1;
            model.nlmsMu = 1;
        }
    }


    import core.runtime : Runtime;
    immutable bool isProgressChecker = Runtime.args.canFind("check");
    immutable bool isResultsMerger = Runtime.args.canFind("merge");


    enum numOfTrials = 10001;
    enum numOfDivTrials = 40;   // N分割する

    static struct SimJob
    {
        ModelSeed modelSeed;
        TaskAppender!(JSONValue, size_t) tasks;
        string dir;
    }

    SimJob[] jobs;

    foreach(dpdMode; [DPDMode.Linearity, DPDMode.EfficiencyAndLinearity])
    foreach(dpdOrder; [1, 7])
    foreach(usePhaseNoise; [false])
    foreach(useIQImbalance; [false])
    foreach(amplifierModel; [/*"Linear",*/ "Rapp", /*"Saleh", /*"Saleh_noPM"*/])
    foreach(numChTaps; [64])
    // ADC&IQ&PA
    foreach(methodName; AliasSeq!(
                                    "L_LS",
                                    // "PHPAOnly3_LS",
                                    // "PHPAOnly5_LS",
                                    "PHPAOnly7_LS",
                                    // "WL_LS",
                                    // "OPH3_LS",
                                    // "OPH5_LS",
                                    // "OPH7_LS",
                                    // "OPHPAOnly7_LS",
                                    // "OPHPAOnly5_LS",
                                    // "OPHPAOnly3_LS",
                                    // "Nop_X",
            ))
    {
        string parentDir = format("results_%sTaps_%s_%s", numChTaps, amplifierModel, numOfTrials);
        if(useIQImbalance)
            parentDir ~= "_withIQI";
        if(usePhaseNoise)
            parentDir ~= "_withPN";
        if(dpdOrder > 1)
            parentDir ~= "_withDPD%s%s".format(cast(string)dpdMode, dpdOrder);

        
        if(dpdOrder == 1 && dpdMode == DPDMode.EfficiencyAndLinearity)
            continue;

        /+
        // only desired signal on AWGN or Rayleigh
        static if(methodName == "Nop_X")
        // foreach(loss_dB; iota(90, 120, 3))
        foreach(txp; iota(10, 32, 1))
        {
            ModelSeed modelSeed;
            modelSeed.txPower = txp.dBm;
            modelSeed.useIQImbalance = useIQImbalance;
            modelSeed.usePhaseNoise = usePhaseNoise;
            modelSeed.useDPD = (dpdOrder > 1);
            modelSeed.orderOfDPD = dpdOrder;
            modelSeed.dpdMode = dpdMode;
            modelSeed.pn3dBBWHz = 0.1;
            modelSeed.paModelName = amplifierModel;
            modelSeed.lnaModelName = amplifierModel;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = 10;
            modelSeed.INR = (-100).dB;
            modelSeed.DesiredLOSS = 100.dB;
            modelSeed.onlyDesired = true;
            modelSeed.outputBER = true;
            modelSeed.linearMode = (amplifierModel == "Linear" ? true : false);
            modelSeed.numOfTapsOfDesiredChannel = numChTaps;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath(parentDir, "results_ber", dir ~ "_onlyDesired");
            dirset[dir] = true;
            appShort.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }
        +/

        // // static if(methodName == "Nop_X" || methodName == "L_LS")
        // foreach(inr_dB; iota(0, 52, 3))
        // {
        //     ModelSeed modelSeed;
        //     modelSeed.linearMode = true;
        //     modelSeed.paModelName = amplifierModel;
        //     modelSeed.lnaModelName = amplifierModel;
        //     modelSeed.paSmoothFactor = paSF;
        //     modelSeed.lnaSmoothFactor = lnaSF;
        //     modelSeed.cancellerType = methodName;
        //     modelSeed.numOfTrainingSymbols = 100;
        //     modelSeed.INR = inr_dB.dB;
        //     modelSeed.SNR = 50.dB;
        //     modelSeed.outputBER = true;
        //     modelSeed.numOfTapsOfDesiredChannel = 64;

        //     auto dir = makeDirNameOfModelSeed(modelSeed);
        //     dir = buildPath(parentDir, "results_ber", dir);
        //     dirset[dir] = true;
        //     appLong.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        // }

        /+
        // INR vs (EVM / SIC / BER)
        foreach(learningSymbols; [200])
        foreach(inr; iota(0, 82, 3))
        foreach(loss; [70])
        {
            ModelSeed modelSeed;
            modelSeed.useIQImbalance = useIQImbalance;
            modelSeed.usePhaseNoise = usePhaseNoise;
            modelSeed.useDPD = (dpdOrder > 1);
            modelSeed.orderOfDPD = dpdOrder;
            modelSeed.dpdMode = dpdMode;
            modelSeed.pn3dBBWHz = 0.1;
            modelSeed.paModelName = amplifierModel;
            modelSeed.lnaModelName = amplifierModel;
            modelSeed.paSmoothFactor = 3;
            modelSeed.lnaSmoothFactor = 3;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = learningSymbols;
            modelSeed.INR = inr.dB;
            // modelSeed.SNR = 50.dB;
            modelSeed.DesiredLOSS = loss.dB;
            modelSeed.onlyDesired = methodName == "Nop_X";
            modelSeed.outputBER = true;
            modelSeed.numofTapsOfSIChannel = numChTaps;
            modelSeed.numOfTapsOfDesiredChannel = numChTaps;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath(parentDir, "results_inr_vs_sic", dir);
            dirset[dir] = true;
            appShort.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }


        // TXP vs (EVM / SIC /  BER)
        foreach(isoRF; [50])
        foreach(desiredLoss; [70/*, 80, 90, 100*/])
        foreach(learningSymbols; [200])
        foreach(txp; iota(10, 32, 1)) {
            ModelSeed modelSeed;
            modelSeed.useIQImbalance = useIQImbalance;
            modelSeed.usePhaseNoise = usePhaseNoise;
            modelSeed.useDPD = (dpdOrder > 1);
            modelSeed.orderOfDPD = dpdOrder;
            modelSeed.dpdMode = dpdMode;
            modelSeed.pn3dBBWHz = 0.1;
            modelSeed.paModelName = amplifierModel;
            modelSeed.lnaModelName = amplifierModel;
            modelSeed.paSmoothFactor = 3;
            modelSeed.lnaSmoothFactor = 3;
            modelSeed.txPower = txp.dBm;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = learningSymbols;
            // modelSeed.INR = (txp-23+inr).dB;
            modelSeed.TXRXISO = isoRF.dB;
            modelSeed.DesiredLOSS = desiredLoss.dB;
            modelSeed.txPower = txp.dBm;
            // modelSeed.onlyDesired = methodName == "Nop_X";
            modelSeed.outputBER = true;
            modelSeed.numofTapsOfSIChannel = numChTaps;
            modelSeed.numOfTapsOfDesiredChannel = numChTaps;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath(parentDir, "results_txp_vs_sic", dir);
            dirset[dir] = true;
            appShort.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }
        +/


        static enum SFLoopConfig
        {
            SAME, X1Only, X2Only
        }


        // SF vs (EVM / SIC / BER)
        if(amplifierModel == "Rapp")
        foreach(txp; [23, 22, 21, 20])
        foreach(isoRF; [50])
        foreach(desiredLoss; [70])
        foreach(learningSymbols; [200])
        foreach(sfconf; [SFLoopConfig.SAME, SFLoopConfig.X1Only, SFLoopConfig.X2Only])
        foreach(sf_10; iota(2, 51, 2)) {
            ModelSeed modelSeed;
            modelSeed.useIQImbalance = useIQImbalance;
            modelSeed.usePhaseNoise = usePhaseNoise;
            modelSeed.useDPD = (dpdOrder > 1);
            modelSeed.orderOfDPD = dpdOrder;
            modelSeed.dpdMode = dpdMode;
            modelSeed.pn3dBBWHz = 0.1;
            modelSeed.paModelName = amplifierModel;
            modelSeed.lnaModelName = amplifierModel;
            // modelSeed.paSmoothFactor = sf_10/10.0;
            modelSeed.paSmoothFactor = 3;
            string saveSuffix;
            final switch(sfconf) {
                case SFLoopConfig.SAME:
                    modelSeed.paSmoothFactor = sf_10/10.0;
                    saveSuffix = "SF1V_SF2V";
                    break;
                case SFLoopConfig.X1Only:
                    modelSeed.paSmoothFactor[0] = sf_10/10.0;
                    saveSuffix = "SF1V_SF2F";
                    break;
                case SFLoopConfig.X2Only:
                    modelSeed.paSmoothFactor[1] = sf_10/10.0;
                    saveSuffix = "SF1F_SF2V";
                    break;
            }

            modelSeed.lnaSmoothFactor = 3;
            modelSeed.txPower = txp.dBm;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = learningSymbols;
            // modelSeed.INR = (txp-23+inr).dB;
            modelSeed.TXRXISO = isoRF.dB;
            modelSeed.DesiredLOSS = desiredLoss.dB;
            modelSeed.onlyDesired = methodName == "Nop_X";
            modelSeed.outputBER = true;
            modelSeed.numofTapsOfSIChannel = numChTaps;
            modelSeed.numOfTapsOfDesiredChannel = numChTaps;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath(parentDir, "results_sf_vs_sic_%s".format(saveSuffix), dir ~ format("_sf%s", sf_10));
            
            auto list = makeTaskListFromModelSeed!methodName(numOfTrials, modelSeed, dir, No.saveAllRAWData);
            jobs ~= SimJob(modelSeed, list, dir);
        }


        // PhaseNoise vs (EVM / SIC / BER)
        // if(usePhaseNoise)
        // foreach(learningSymbols; [200])
        // foreach(pn3dBBWHz; iota(-20, 23, 2).map!(a => 10.0^^(a/10.0)))
        // foreach(loss; [70])
        // {
        //     ModelSeed modelSeed;
        //     modelSeed.pn3dBBWHz = pn3dBBWHz;
        //     modelSeed.useIQImbalance = useIQImbalance;
        //     modelSeed.usePhaseNoise = usePhaseNoise;
        //     modelSeed.useDPD = (dpdOrder > 1);
        //     modelSeed.orderOfDPD = dpdOrder;
            // modelSeed.dpdMode = dpdMode;
        //     modelSeed.paModelName = amplifierModel;
        //     modelSeed.lnaModelName = amplifierModel;
        //     modelSeed.paSmoothFactor = 3;
        //     modelSeed.lnaSmoothFactor = 3;
        //     modelSeed.cancellerType = methodName;
        //     modelSeed.numOfTrainingSymbols = learningSymbols;
        //     modelSeed.INR = 60.dB;
        //     // modelSeed.SNR = 50.dB;
        //     modelSeed.DesiredLOSS = loss.dB;
            // modelSeed.onlyDesired = methodName == "Nop_X";
        //     modelSeed.outputBER = true;
        //     modelSeed.numofTapsOfSIChannel = numChTaps;
        //     modelSeed.numOfTapsOfDesiredChannel = numChTaps;

        //     auto dir = makeDirNameOfModelSeed(modelSeed);
        //     dir = buildPath(parentDir, "results_phasenoise_vs_sic", dir);
        //     dirset[dir] = true;
        //     appShort.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        // }
    }


    SplitMergeResumableTasks!(typeof(SimJob.init.tasks))[] splittedTaskList;
    foreach(i, jb; jobs) {
        if(thisProcessType() == ChildProcessType.SUBMITTER && !isProgressChecker && !isResultsMerger)
            mkdirRecurse(jb.dir);

        auto sp = jb.tasks.toSplitMergeResumable(numOfDivTrials, buildPath(jb.dir, "savedata"), 50);
        splittedTaskList ~= sp;
    }

    import std.stdio;
    import core.runtime : Runtime;
    import std.digest.crc : crc32Of;
    import std.digest : toHexString;
    import std.datetime;

    enforce(!isProgressChecker);

    if(isResultsMerger)
    {
        JobEnvironment env = defaultJobEnvironment();
        env.maxArraySize = 1;
        env.jobName = format("merge_%s_%s", Runtime.args[0].baseName, crc32Of(cast(ubyte[])std.file.read(Runtime.args[0])).toHexString);

        auto list = new MultiTaskList!void();
        foreach(i, sp; splittedTaskList) {
            list.append((){
                auto dir = jobs[i].dir;
                auto resList = sp.returns;
                sp.nullify();
                enforce(resList.length == numOfTrials, "Tasks (%s) have not yet been completed.".format(dir));

                JSONValue jv = cast(JSONValue[string])null;
                jv["RemainPowers"] = resList.map!(a => a["RemainPower"]).array();
                jv["SIPowers"] = resList.map!(a => a["SIPower"]).array();
                jv["cancellations"] = resList.map!(a => a["cancellation_dB"]).array();
                jv["RINRs"] = resList.map!(a => a["RINR_dB"]).array();
                jv["INRs"] = resList.map!(a => a["INR_dB"]).array();
                jv["TXPs"] = resList.map!(a => a["TXPower_dBm"]).array();
                if(jobs[i].modelSeed.outputBER){
                    jv["bers"] = resList.map!(a => a["ber"]).array();
                    // jv["evms"] = resList.map!(a => a["evm"]).array();
                    jv["sers"] = resList.map!(a => a["ser"]).array();
                }

                writeJSONData(buildPath(dir, ALL_RESULT_FILENAME), jv);
            });
        }

        tuthpc.taskqueue.run(list, env);
    }
    else
    {
        auto taskList = new MultiTaskList!void();
        foreach(sp; splittedTaskList)
            taskList ~= sp.toMultiTaskList();

        JobEnvironment env = defaultJobEnvironment();
        env.pmem = 2;
        env.mem = env.pmem * env.taskGroupSize;
        env.jobName = format("run_%s_%s", Runtime.args[0].baseName, crc32Of(cast(ubyte[])std.file.read(Runtime.args[0])).toHexString);

        if(taskList.length != 0)
            tuthpc.taskqueue.run(taskList, env);
    }
}


Model[] makeModels(string methodName)(size_t numOfTrials, ModelSeed modelSeed, string uniqueString, size_t startIndex = 0, size_t endIndex = size_t.max)
{
    import dffdd.blockdiagram.noise : noisePower;
    Model[] models;

    immutable string hashOfModelSeed = () @trusted {
        import std.digest.crc;
        import tuthpc.taskqueue : hashOfExe;
        return crc32Of(methodName, (cast(ubyte*)&modelSeed)[0 .. modelSeed.sizeof], hashOfExe(), uniqueString).crcHexString;
    }();

    foreach(iTrial; startIndex .. min(numOfTrials, endIndex)) {
        Model model;
        scope(success)
            models ~= model;

        immutable NP = (10*log10(noisePower(model.samplingFreq, 300)) + 30).dBm * 4.dB;

        model.uniqueId = format("%s_%s", hashOfModelSeed, iTrial.to!string);
        model.rndSeed = (cast(uint)hashOf(iTrial));
        // model.SNR = modelSeed.SNR;

        if(modelSeed.SNR.isNull)
            model.SNR = modelSeed.txPower[0] / NP / modelSeed.DesiredLOSS.get;
        else
            model.SNR = modelSeed.SNR.get;

        if(modelSeed.INR.isNull)
            model.INR = modelSeed.txPower[0] / NP / modelSeed.TXRXISO.get;
        else
            model.INR = modelSeed.INR.get;

        // model.pa.TX_POWER = modelSeed.txPower;
        // model.txIQMixer.IRR = modelSeed.txIRR;
        // model.rxIQMixer.IRR = modelSeed.rxIRR;
        model.basisFuncsSelection.noiseMargin = modelSeed.gamma;
        foreach(i; 0 .. 2) model.xcvrs[i].lna.smoothFactor = modelSeed.lnaSmoothFactor[i];

        model.withSI = !modelSeed.onlyDesired;

        // model.withSI = false;
        // model.withSIC = false;

        // 再現する非線形性の選択
        model.useDTXDPD = modelSeed.useDPD;
        model.useSTXDPD = modelSeed.useDPD;
        model.useDTXIQ = modelSeed.useIQImbalance;
        // model.useDTXPN = true;
        model.useDTXPN = modelSeed.usePhaseNoise;
        model.useDTXPA = ! modelSeed.linearMode;
        model.useSTXIQ = modelSeed.useIQImbalance;
        model.useSTXPN = modelSeed.usePhaseNoise;
        model.useSTXPA = ! modelSeed.linearMode;
        model.useSRXLN = ! modelSeed.linearMode;
        model.useSRXIQ = modelSeed.useIQImbalance;
        model.useSRXQZ = false;


        /* PAの設定 */
        foreach(i; 0 .. 2)
        {
            model.xcvrs[i].pa.modelName = modelSeed.paModelName[i];
            if(model.xcvrs[i].pa.modelName == "Linear") {
                model.useSTXPA = false;
                model.useDTXPA = false;
            }

            model.xcvrs[i].pa.TX_POWER = modelSeed.txPower[i];
            model.xcvrs[i].pa.GAIN = modelSeed.paGain[i];
            model.xcvrs[i].pa.Vsat = modelSeed.paVsat[i];
            model.xcvrs[i].pa.smoothFactor = modelSeed.paSmoothFactor[i];
        }

        /* LNAの設定 */
        foreach(i; 0 .. 2)
        {
            model.xcvrs[i].lna.modelName = modelSeed.lnaModelName[i];
            if(modelSeed.lnaModelName[i] == "Linear") {
                model.useSRXLN = false;
            }

            model.xcvrs[i].lna.GAIN = modelSeed.lnaGain[i];
            model.xcvrs[i].lna.Vsat = modelSeed.lnaVsatIn[i] * modelSeed.lnaGain[i];   // 入力飽和電圧から出力飽和電圧への変換
            model.xcvrs[i].lna.smoothFactor = modelSeed.lnaSmoothFactor[i];
        }

        // 位相雑音
        model.phaseNoise.betaBWHz = modelSeed.pn3dBBWHz;

        /* TX IQ Mixer の設定 */
        foreach(i; 0 .. 2)
        {
            Random rnd = uniqueRandom(iTrial, "TXIQMixer_%s".format(i));
            model.xcvrs[i].txIQMixer.imbCoef = (1.0L / modelSeed.txIRR.asV) * std.complex.expi(uniform(0, 1.0f, rnd) * 2 * PI);
        }

        /* RX IQ Mixer の設定 */
        foreach(i; 0 .. 2)
        {
            Random rnd = uniqueRandom(iTrial, "RXIQMixer_%s".format(i));
            model.xcvrs[i].rxIQMixer.imbCoef = (1.0L / modelSeed.rxIRR.asV) * std.complex.expi(uniform(0, 1.0f, rnd) * 2 * PI);
        }

        /* チャネルの設定 */
        {
            Random rnd = uniqueRandom(iTrial, "Channel");
            model.channelSI.taps = modelSeed.numofTapsOfSIChannel;

            BoxMuller!Random gGen = BoxMuller!Random(rnd);
            Complex!double[] coefs;
            foreach(i; 0 .. model.channelSI.taps){
                // tapsのタップ数で40dB減衰
                // auto db = -1 * (40.0 / model.channelSI.taps) * i;
                auto db = 0;
                coefs ~= cast(Complex!double)(gGen.front * 10.0 ^^ (db/20));
                gGen.popFront();
            }

            model.channelSI.impulseResponse = coefs;
        }
        {
            Random rnd = uniqueRandom(iTrial, "ChannelDesired");
            model.channelDesired.taps = modelSeed.numOfTapsOfDesiredChannel;

            BoxMuller!Random gGen = BoxMuller!Random(rnd);
            Complex!double[] coefs;
            foreach(i; 0 .. model.channelDesired.taps) {
                // tapsのタップ数で40dB減衰
                // auto db = -1 * (40 / (model.ofdm.numOfCP * model.ofdm.scaleOfUpSampling)) * i;
                auto db = 0;
                coefs ~= cast(Complex!double)(gGen.front * 10.0 ^^ (db/20));
                gGen.popFront();
            }

            model.channelDesired.impulseResponse = coefs;
        }

        /* キャンセラの設定 */
        {
            model.orthogonalizer.numOfTrainingSamples = 10_000;
            model.firFilter.taps = model.channelSI.taps;

            if(methodName[0] == 'O')
                model.orthogonalizer.enabled = true;
            else
                model.orthogonalizer.enabled = false;

            // ベースバンド信号波形の出力
            model.outputWaveform = false;

            // BERの出力
            model.outputBER = modelSeed.outputBER;

            if(methodName.canFind("DCM")) {
                model.learningSymbols = 10;
                model.learningCount = 3;
            } else {
                model.learningSymbols = modelSeed.numOfTrainingSymbols;
                model.learningCount = 1;
            }

            model.numOfFilterTrainingSymbols = max(model.learningSymbols * model.learningCount * 2 + 10, 100);

            model.swappedSymbols = 100000;
            model.rlsAdapter.delta = modelSeed.rlsDelta;
            model.rlsAdapter.lambda = modelSeed.rlsLambda;
            model.nlmsAdapter.mu = modelSeed.nlmsMu;
            model.iterativeFreqSIC.iterations = modelSeed.iterNumOfIteration;
            model.iterativeFreqSIC.newtonIterations = modelSeed.iterNumOfNewton;
            model.iterativeFreqSIC.numOfSCForEstNL = modelSeed.iterNumOfSCForEstNL;
            model.iterativeFreqSIC.estimationOrder = modelSeed.iterEstOrder;
            model.iterativeFreqSIC.use3rdSidelobe = modelSeed.iterUse3rdSidelobe;
        }


        /* DPDの設定 */
        {
            model.dpd.order = modelSeed.orderOfDPD;
            model.dpd.mode = modelSeed.dpdMode;
            model.dpd.numOfTrainingSymbols = 10;
        }
    }

    return models;
}


string makeDirNameOfModelSeed(ModelSeed modelSeed)
{
    string strINRorISO;
    if(modelSeed.INR.isNull)
        strINRorISO = format("iso%s", modelSeed.TXRXISO.get);
    else
        strINRorISO = format("inr%s", modelSeed.INR.get);

    string strSNRorISO;
    if(modelSeed.SNR.isNull)
        strSNRorISO = format("loss%s", modelSeed.DesiredLOSS.get);
    else
        strSNRorISO = format("snr%s", modelSeed.SNR.get);

    string cancellerName = modelSeed.cancellerType;
    if(cancellerName.canFind("Sidelobe"))
        cancellerName ~= "_" ~ modelSeed.iterEstOrder;

    string dir;
    if(modelSeed.cancellerType.split("_")[0].endsWith("FHF"))
    {
        dir = "TXP%s_%s_%s_%s_%s_G%s_IRR%s_%s"
            .format(modelSeed.txPower[0], modelSeed.txPower[1],
                    strINRorISO,
                    strSNRorISO,
                    cancellerName,
                    modelSeed.bfsNoiseMargin,
                    modelSeed.txIRR,
                    modelSeed.numOfTrainingSymbols);
    }
    else
    {
        dir = "TXP%s_%s_%s_%s_%s_IRR%s_PN%s_%s"
            .format(modelSeed.txPower[0], modelSeed.txPower[1],
                    strINRorISO,
                    strSNRorISO,
                    cancellerName,
                    modelSeed.txIRR,
                    10*log10(modelSeed.pn3dBBWHz)+100,
                    modelSeed.numOfTrainingSymbols);
    }

    return dir;
}


Random uniqueRandom(Args...)(Args args)
{
    ulong v = 0;
    foreach(ref p; args)
        v += hashOf(p);

    Random rnd;
    rnd.seed(v & uint.max);
    return rnd;
}


JSONValue runTrial(string methodName)(size_t indexTrial, Model m, string dir, Flag!"saveAllRAWData" saveAllRAWData = No.saveAllRAWData)
{
    import core.memory;
    GC.collect();
    GC.minimize();

    JSONValue res;
    if(indexTrial == 0 && !NO_NEED_CSV_FILES)
        res = mainImpl!(methodName)(m,  dir);
    else
        res = mainImpl!(methodName)(m,  null);

    res["model"] = (){
        static JSONValue cpx2JV(F)(Complex!F a) { return JSONValue(["re": a.re, "im": a.im]); }

        JSONValue jv = cast(JSONValue[string])null;
        jv["impulseResponseSI"] = m.channelSI.impulseResponse.map!cpx2JV.array();
        jv["impulseResponseDesired"] = m.channelDesired.impulseResponse.map!cpx2JV.array();

        foreach(i; 0 .. 2) {
            jv["txIQCoef_%s".format(i)] = cpx2JV(m.xcvrs[i].txIQMixer.imbCoef);
            jv["rxIQCoef_%s".format(i)] = cpx2JV(m.xcvrs[i].rxIQMixer.imbCoef);
        }

        return jv;
    }();

    return res;
}


TaskAppender!(JSONValue, size_t) makeTaskListFromModelSeed(string methodName)(size_t nTrials, ModelSeed modelSeed, string dir, Flag!"saveAllRAWData" saveAllRAWData = No.saveAllRAWData)
{
    JSONValue runImpl(size_t i)
    {
        auto model = makeModels!methodName(nTrials, modelSeed, dir, i, i+1)[0];
        return runTrial!methodName(i, model, dir, saveAllRAWData);
    }

    auto app = taskAppender(&runImpl);
    foreach(i; 0 .. nTrials) {
        app.append(i);
    }

    return app;
}


void main()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    // jobRun(1, 0, {
        mainJob();
    // });
}
