module app;

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

import models;
import simmain;

import msgpack;

extern(C) void openblas_set_num_threads(int num_threads);
extern(C) int openblas_get_num_threads();


enum NO_NEED_CSV_FILES = true;
immutable string ALL_RESULT_FILENAME = "all_result.bin";
immutable string DUMPED_RESULT_LIST_FILENAME = "dumped_results.bin";

void writeJSONData(string path, JSONValue jv)
{
    std.file.write(path, msgpack.pack(msgpack.fromJSONValue(jv)));
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

    /* IQ Mixer */
    Gain txIRR = 25.dB;
    Gain rxIRR = 25.dB;

    /* PA */
    // Voltage txPower = 23.dBm;
    // Voltage paIIP3 = 20.dBm;
    // Gain paGain = 27.dB;
    enum real txBackoff_dB = 7;
    enum real txPower_dBm = 23;
    enum real paGain_dB = 30;
    Voltage txPower = txPower_dBm.dBm;
    Voltage paIIP3 = (txPower_dBm - paGain_dB + 6 + txBackoff_dB).dBm;
    Gain paGain = paGain_dB.dB;
    double paSmoothFactor = 3;

    /* LNA */
    double lnaSmoothFactor = 3;

    /* Other */
    Nullable!Gain SNR;
    Nullable!Gain DesiredLOSS;
    Nullable!Gain INR;
    Nullable!Gain TXRXISO;
    uint numOfTrainingSymbols = 20;
    Gain gamma = 0.dB;

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

    auto taskListShort = new MultiTaskList!void();
    auto taskListLong = new MultiTaskList!void();

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


    enum bool isProgressChecker = true;
    enum bool isDumpedFileCheck = true;


    enum numOfTrials = 101;
    size_t sumOfTaskNums = 0;
    size_t sumOfTrials = 0;


    foreach(numChTaps; [64])
    // ADC&IQ&PA
    foreach(methodName; AliasSeq!(
                                    "OPHPAOnly3_LS",
                                    "OPHPAOnly5_LS",
                                    "OPHPAOnly7_LS",
                                    // "OPHPAOnly5_LS",
                                    // "OPHPAOnly3_LS",
                                    // "Nop_X",

                                    "L_LS",
            ))
    {
        bool[string] dirset;
        auto appShort = uniqueTaskAppender(&mainForEachTrial!methodName);
        auto appLong = uniqueTaskAppender(&mainForEachTrial!methodName);
        scope(exit){
            enforce(appShort.length + appLong.length == dirset.length);
            taskListShort ~= appShort;
            taskListLong ~= appLong;
            
            static if(isProgressChecker)
            {
                foreach(dir, _; dirset) {
                    if(!exists(buildPath(dir, ALL_RESULT_FILENAME))){
                        sumOfTaskNums += 1;
                        // writeln(dir);

                        sumOfTrials += numOfTrials;

                        string resListDumpFilePath = buildPath(dir, DUMPED_RESULT_LIST_FILENAME);
                        if(isDumpedFileCheck && exists(resListDumpFilePath)) {
                            // size_t compls = std.file.readText(resListDumpFilePath)
                            //                 .parseJSON(JSONOptions.specialFloatLiterals).array.length;
                            // size_t compls = msgpack.toJSONValue(msgpack.unpack(cast(ubyte[]) std.file.read(resListDumpFilePath))).array.length;
                            size_t compls = readJSONData(resListDumpFilePath).array.length;
                            sumOfTrials -= compls;
                            writefln!"%s: %s complete"(dir, compls);
                        }
                    }
                }
            }
        }


        string parentDir = format("results_Taps%s", numChTaps);


        // only desired signal on AWGN or Rayleigh
        static if(methodName == "Nop_X")
        foreach(snr; iota(0, (numChTaps == 1 ? 22 : 52), 3))
        {
            ModelSeed modelSeed;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = 10;
            modelSeed.INR = 20.dB;
            modelSeed.SNR = snr.dB;
            modelSeed.onlyDesired = true;
            modelSeed.outputBER = true;
            modelSeed.linearMode = true;
            modelSeed.numOfTapsOfDesiredChannel = numChTaps;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath(parentDir, "results_ber", dir ~ "_onlyDesired");
            dirset[dir] = true;
            appShort.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }

        // // static if(methodName == "Nop_X" || methodName == "L_LS")
        // foreach(inr_dB; iota(0, 52, 3))
        // {
        //     ModelSeed modelSeed;
        //     // modelSeed.linearMode = true;
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


        // INR vs (EVM / SIC / BER)
        foreach(learningSymbols; [50])
        foreach(inr; iota(0, 82, 3))
        foreach(loss; [70])
        {
            ModelSeed modelSeed;
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
        foreach(isoRF; [50, 70])
        foreach(desiredLoss; [70])
        foreach(learningSymbols; [50])
        foreach(txp; iota(10, 32, 1)) {
            ModelSeed modelSeed;
            modelSeed.paSmoothFactor = 3;
            modelSeed.lnaSmoothFactor = 3;
            modelSeed.txPower = txp.dBm;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = learningSymbols;
            // modelSeed.INR = (txp-23+inr).dB;
            modelSeed.TXRXISO = isoRF.dB;
            modelSeed.DesiredLOSS = desiredLoss.dB;
            modelSeed.txPower = txp.dBm;
            modelSeed.onlyDesired = methodName == "Nop_X";
            modelSeed.outputBER = true;
            modelSeed.numofTapsOfSIChannel = numChTaps;
            modelSeed.numOfTapsOfDesiredChannel = numChTaps;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath(parentDir, "results_txp_vs_sic", dir);
            dirset[dir] = true;
            appShort.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }


        // SF vs (EVM / SIC / BER)
        foreach(txp; [23])
        foreach(isoRF; [50, 70])
        foreach(desiredLoss; [70])
        foreach(learningSymbols; [50])
        foreach(sf_10; iota(2, 51, 2)) {
            ModelSeed modelSeed;
            modelSeed.paSmoothFactor = sf_10/10.0;
            modelSeed.lnaSmoothFactor = sf_10/10.0;
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
            dir = buildPath(parentDir, "results_sf_vs_sic", dir ~ format("_sf%s", sf_10));
            dirset[dir] = true;
            appShort.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }
    }

    import std.stdio;


    static if(isProgressChecker)
    {
        immutable size_t totalTaskListLen = taskListShort.length + taskListLong.length;

        writefln("%s tasks will be submitted.", totalTaskListLen);
        writefln("%s tasks will be computed.", sumOfTaskNums);

        immutable size_t totalTrials = numOfTrials * totalTaskListLen;
        immutable size_t completeTrials = totalTrials - sumOfTrials;
        writefln("%s/%s (%2.2f%%) is completed. ", completeTrials, totalTrials, completeTrials*1.0/totalTrials*100);
    }
    else
    {
        JobEnvironment env = defaultJobEnvironment();
        env.pmem = 3;
        env.mem = env.pmem * env.taskGroupSize;
        // env.scriptPath = "jobscript.sh";

        if(taskListShort.length != 0)
            tuthpc.taskqueue.run(taskListShort, env);

        if(taskListLong.length != 0)
            tuthpc.taskqueue.run(taskListLong, env);
    }

    // foreach(i; 0 .. taskListShort.length)
    //     taskListShort[i]();

    // foreach(i; 0 .. taskListLong.length)
    //     taskListLong[i]();
}


Model[] makeModels(string methodName)(size_t numOfTrials, ModelSeed modelSeed, string uniqueString)
{
    import dffdd.blockdiagram.noise : noisePower;
    Model[] models;

    immutable string hashOfModelSeed = () @trusted {
        import std.digest.crc;
        import tuthpc.taskqueue : hashOfExe;
        return crc32Of(methodName, (cast(ubyte*)&modelSeed)[0 .. modelSeed.sizeof], hashOfExe(), uniqueString).crcHexString;
    }();

    foreach(iTrial; 0 .. numOfTrials) {
        Model model;
        scope(success)
            models ~= model;

        immutable NP = (10*log10(noisePower(model.samplingFreq, 300)) + 30).dBm * 4.dB;

        model.uniqueId = format("%s_%s", hashOfModelSeed, iTrial.to!string);
        model.rndSeed = cast(uint)hashOf(iTrial);
        // model.SNR = modelSeed.SNR;

        if(modelSeed.SNR.isNull)
            model.SNR = modelSeed.txPower / NP / modelSeed.DesiredLOSS;
        else
            model.SNR = modelSeed.SNR;

        if(modelSeed.INR.isNull)
            model.INR = modelSeed.txPower / NP / modelSeed.TXRXISO;
        else
            model.INR = modelSeed.INR;

        // model.pa.TX_POWER = modelSeed.txPower;
        // model.txIQMixer.IRR = modelSeed.txIRR;
        // model.rxIQMixer.IRR = modelSeed.rxIRR;
        model.basisFuncsSelection.noiseMargin = modelSeed.gamma;
        model.lna.smoothFactor = modelSeed.lnaSmoothFactor;

        model.withSI = !modelSeed.onlyDesired;

        // model.withSI = false;
        // model.withSIC = false;

        // 再現する非線形性の選択
        model.useDTXIQ = false;
        model.useDTXPN = false;
        model.useDTXPA = ! modelSeed.linearMode;
        model.useSTXIQ = false;
        model.useSTXPN = false;
        model.useSTXPA = ! modelSeed.linearMode;
        model.useSRXLN = ! modelSeed.linearMode;
        model.useSRXIQ = false;
        model.useSRXQZ = false;


        /* PAの設定 */
        {
            model.pa.TX_POWER = modelSeed.txPower;
            model.pa.IIP3 = modelSeed.paIIP3;
            model.pa.GAIN = modelSeed.paGain;
            model.pa.smoothFactor = modelSeed.paSmoothFactor;
        }

        /* TX IQ Mixer の設定 */
        {
            Random rnd = uniqueRandom(iTrial, "TXIQMixer");
            model.txIQMixer.imbCoef = model.txIQMixer.imbCoef = (1.0L / modelSeed.txIRR.asV) * std.complex.expi(uniform(0, 1.0f, rnd) * 2 * PI);
        }

        /* RX IQ Mixer の設定 */
        {
            Random rnd = uniqueRandom(iTrial, "RXIQMixer");
            model.rxIQMixer.imbCoef = model.txIQMixer.imbCoef = (1.0L / modelSeed.rxIRR.asV) * std.complex.expi(uniform(0, 1.0f, rnd) * 2 * PI);
        }

        /* チャネルの設定 */
        {
            Random rnd = uniqueRandom(iTrial, "Channel");
            model.channelSI.taps = 64;

            BoxMuller!Random gGen = BoxMuller!Random(rnd);
            Complex!real[] coefs;
            foreach(i; 0 .. model.channelSI.taps){
                // tapsのタップ数で40dB減衰
                auto db = -1 * (40.0 / model.channelSI.taps) * i;
                coefs ~= cast(Complex!real)(gGen.front * 10.0L ^^ (db/20));
                gGen.popFront();
            }

            model.channelSI.impulseResponse = coefs;
        }
        {
            Random rnd = uniqueRandom(iTrial, "ChannelDesired");
            model.channelDesired.taps = modelSeed.numOfTapsOfDesiredChannel;

            BoxMuller!Random gGen = BoxMuller!Random(rnd);
            Complex!real[] coefs;
            foreach(i; 0 .. model.channelDesired.taps) {
                // tapsのタップ数で40dB減衰
                auto db = -1 * (40 / (model.ofdm.numOfCP * model.ofdm.scaleOfUpSampling)) * i;
                coefs ~= cast(Complex!real)(gGen.front * 10.0L ^^ (db/20));
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
        dir = "TXP%s_%s_%s_%s_G%s_IRR%s_%s"
            .format(modelSeed.txPower,
                    strINRorISO,
                    strSNRorISO,
                    cancellerName,
                    modelSeed.bfsNoiseMargin,
                    modelSeed.txIRR,
                    modelSeed.numOfTrainingSymbols);
    }
    else
    {
        dir = "TXP%s_%s_%s_%s_IRR%s_%s"
            .format(modelSeed.txPower,
                    strINRorISO,
                    strSNRorISO,
                    cancellerName,
                    modelSeed.txIRR,
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


void mainForEachTrial(string methodName)(size_t nTrials, ModelSeed modelSeed, string dir, Flag!"saveAllRAWData" saveAllRAWData = No.saveAllRAWData)
{
    Model[] models = makeModels!methodName(nTrials, modelSeed, dir);

    if(exists(buildPath(dir, ALL_RESULT_FILENAME))) return;

    writeln(dir);
    mkdirRecurse(dir);

    // 前回異常終了等で死んでいれば，このファイルがあるはず
    // 中間状態のダンプファイル
    immutable resListDumpFilePath = buildPath(dir, DUMPED_RESULT_LIST_FILENAME);
    JSONValue[] lastDumpedList;
    if(exists(resListDumpFilePath))
    {
        // lastDumpedList = std.file.readText(resListDumpFilePath)
        //     .parseJSON(JSONOptions.specialFloatLiterals).array;
        // lastDumpedList = msgpack.toJSONValue(msgpack.unpack(cast(ubyte[]) std.file.read(resListDumpFilePath))).array;
        lastDumpedList = readJSONData(resListDumpFilePath).array;
    }
    scope(success) {
        // このプロセスが正常終了すれば中間状態のダンプファイルは不要
        if(exists(resListDumpFilePath))
            std.file.remove(buildPath(dir, DUMPED_RESULT_LIST_FILENAME));
    }

    import std.datetime;
    auto lastUpdateTime = Clock.currTime;


    JSONValue[] resList;
    uint sumOfSuccFreq;
    JSONValue[] selectingRatioList;
    foreach(i, ref m; models) {
        {
            import core.memory;
            GC.collect();
        }

        JSONValue res; // この試行での結果が格納される
        scope(success)
        {
            resList ~= res;

            // このプロセスで計算できた試行回数が前回を上回っていればダンプファイルを更新
            // ただし，前回の更新時より1分以上空いていること
            auto ct = Clock.currTime;
            if(resList.length > lastDumpedList.length
                && (ct - lastUpdateTime) > 60.seconds )
            {
                // std.file.write(resListDumpFilePath, JSONValue(resList).toString(JSONOptions.specialFloatLiterals));
                // File outputFile = File(resListDumpFilePath, "w");
                // auto packer = msgpack.Packer(outputFile.lockingTextWriter, false);
                // packer.pack(msgpack.fromJSONValue(JSONValue(resList)));
                // std.file.write(resListDumpFilePath, msgpack.pack(msgpack.fromJSONValue(JSONValue(resList))));
                writeJSONData(resListDumpFilePath, JSONValue(resList));
                lastUpdateTime = ct;
            }
        }

        if(lastDumpedList.length >= i+1) {
            // 前回中断時にすでに計算済みなのでそのデータを復元
            res = lastDumpedList[i];
            continue;
        }
        else {
            // まだ未計算なので計算する
            if(i == 0 && !NO_NEED_CSV_FILES)
                res = mainImpl!(methodName)(m,  dir);
            else
                res = mainImpl!(methodName)(m,  null);
        }

        if(saveAllRAWData) {
            res["model"] = (){
                static JSONValue cpx2JV(F)(Complex!F a) { return JSONValue(["re": a.re, "im": a.im]); }

                JSONValue jv = cast(JSONValue[string])null;
                jv["impulseResponseSI"] = m.channelSI.impulseResponse.map!cpx2JV.array();
                jv["impulseResponseDesired"] = m.channelDesired.impulseResponse.map!cpx2JV.array();
                jv["txIQCoef"] = cpx2JV(m.txIQMixer.imbCoef);
                jv["rxIQCoef"] = cpx2JV(m.rxIQMixer.imbCoef);

                return jv;
            }();
        }
    }

    JSONValue jv = cast(JSONValue[string])null;
    jv["cancellations"] = resList.map!(a => a["cancellation_dB"]).array();
    jv["RINRs"] = resList.map!(a => a["RINR_dB"]).array();
    jv["INRs"] = resList.map!(a => a["INR_dB"]).array();
    if(modelSeed.outputBER){
        jv["bers"] = resList.map!(a => a["ber"]).array();
        jv["evms"] = resList.map!(a => a["evm"]).array();
        jv["sers"] = resList.map!(a => a["ser"]).array();
    }

    // auto file = File(buildPath(dir, "allResult.json"), "w");
    // file.write(jv.toPrettyString(JSONOptions.specialFloatLiterals));
    {
        // auto file = File(buildPath(dir, ALL_RESULT_FILENAME), "w");
        // file.rawWrite(msgpack.pack(msgpack.fromJSONValue(jv)));
        writeJSONData(buildPath(dir, ALL_RESULT_FILENAME), jv);
    }

    if(saveAllRAWData){
        // auto file = File(buildPath(dir, "rawAllResult.bin"), "w");
        // file.write(JSONValue(resList).toPrettyString(JSONOptions.specialFloatLiterals));
        // file.rawWrite(msgpack.pack(msgpack.fromJSONValue(JSONValue(resList))));
        writeJSONData(buildPath(dir, "rawAllResult.bin"), JSONValue(resList));
    }

    if(methodName.startsWith("Log_")){
        foreach(i, res; resList){
            import std.file : mkdirRecurse;

            auto siglogpath = buildPath(dir, "signallog");
            mkdirRecurse(siglogpath);

            auto uniqueId = res["filterSpec"]["uniqueId"].str;
            foreach(iden; ["TrX", "TrY", "TrZ", "CnX", "CnY", "CnZ"]){
                string filename = "signal_%s_%s.dat".format(iden, uniqueId);
                auto frompath = buildPath("signallog", filename);
                auto topath = buildPath(siglogpath, "signal_%s_%s.dat".format(iden, i));
                
                // frompathからtopathにファイルをコピーする
                std.file.rename(frompath, topath);
                writefln("%s -> %s", frompath, topath);
            }
        }
    }
}


void main()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    // jobRun(1, 0, {
        mainJob();
    // });
}
