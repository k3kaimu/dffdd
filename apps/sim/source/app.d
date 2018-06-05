module app;

import std.algorithm;
import std.format;
import std.json;
import std.math;
import std.meta;
import std.path;
import std.range;
import std.stdio;
import std.random;
import std.range;
import std.file;
import std.typecons;
import std.exception;
import std.complex;

import dffdd.utils.unit;
import dffdd.blockdiagram.noise;

import models;
import simmain;


extern(C) void openblas_set_num_threads(int num_threads);
extern(C) int openblas_get_num_threads();


shared static this()
{
    openblas_set_num_threads(1);
}



struct ModelSeed
{
    string cancellerType;

    /* IQ Mixer */
    Gain txIRR = 25.dB;
    Gain rxIRR = 25.dB;

    /* PA */
    Voltage txPower = 23.dBm;
    Voltage paIIP3 = 21.8.dBm;
    Gain paGain = 28.5.dB;

    /* LNA */
    uint lnaSmoothFactor = 3;

    /* Other */
    Gain SNR = 11.dB;
    Gain INR = 50.dB;
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
    uint iterNumOfIteration = 10;

    /* measure only desired signal */
    bool onlyDesired = false;
}



void mainJob()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    auto taskList = new MultiTaskList();

    auto setRLSLMSParam(ref ModelSeed model)
    {
        if(model.cancellerType.startsWith("OPH")) {
            model.rlsDelta = 0.1;
            model.rlsLambda = 1 - 1.5E-3;
            model.nlmsMu = 0.33;
        } else {
            model.rlsDelta = 4E-7;
            model.rlsLambda = 1;
            model.nlmsMu = 1;
        }
    }


    enum numOfTrials = 1;

    // ADC&IQ&PA
    foreach(methodName; AliasSeq!(
                                    // "S2FHF_LS",
                                    // "S2FHF_RLS",
                                    // "S2FHF_LMS",
                                    // //
                                    // "FHF_LS",
                                    // "FHF_RLS",
                                    // "FHF_LMS",
                                    // //
                                    // "OPH_RLS",
                                    "OPH_LS",
                                    // "OPH_LMS",
                                    
                                    // "WL_LS",
                                    // "L_LS",
                                    
                                    // "IterativeFreqSIC_X",
                                    // "SidelobeFwd_X",
                                    // "SidelobeInv_X"
            ))
    {
        bool[string] dirset;
        auto appender = uniqueTaskAppender(&mainForEachTrial!methodName);
        scope(exit){
            enforce(appender.length == dirset.length);
            taskList ~= appender;
        }


        /// inr vs cancellation
        // static if(methodName.endsWith("_LS") || methodName == "IterativeFreqSIC_X")
        foreach(learningSymbols; [100])
        // foreach(inr; iota(0, 12, 2))
        foreach(snr; iota(0, 21, 5))
        foreach(sf; [1])
        foreach(siflag; [false, true])
        {
            ModelSeed modelSeed;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = learningSymbols;
            modelSeed.INR = 20.dB;
            modelSeed.SNR = snr.dB;
            modelSeed.onlyDesired = siflag;
            import dffdd.mod.qam;
            writefln("BER = %s at SNR = %sdB", berQAMFromSNR(snr, 16), snr);

            modelSeed.lnaSmoothFactor = sf;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath("results_ALL_inr_vs_canc", dir ~ "_%s".format(siflag ? "withSI" : "onlyDesired"));
            dirset[dir] = true;

            appender.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }
    }

    import std.stdio;

    //writefln("%s tasks will be submitted.", taskList.length);
    JobEnvironment env;
    tuthpc.taskqueue.run(taskList, env);
    // foreach(i; 0 .. taskList.length)
        // taskList[i]();
}


Model[] makeModels(string methodName)(size_t numOfTrials, ModelSeed modelSeed)
{
    Model[] models;

    foreach(iTrial; 0 .. numOfTrials) {
        Model model;
        scope(success)
            models ~= model;

        model.rndSeed = cast(uint)hashOf(iTrial);
        model.SNR = modelSeed.SNR;
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
        model.useDTXPA = false;
        model.useSTXIQ = true;
        model.useSTXPN = false;
        model.useSTXPA = true;
        model.useSRXLN = true;
        model.useSRXIQ = true;
        model.useSRXQZ = true;


        /* PAの設定 */
        {
            model.pa.TX_POWER = modelSeed.txPower;
            model.pa.IIP3 = modelSeed.paIIP3;
            model.pa.GAIN = modelSeed.paGain;
        }

        /* TX IQ Mixer の設定 */
        {
            Random rnd = uniqueRandom(iTrial, "TXIQMixer");
            model.txIQMixer.imbCoef = (1.0L / modelSeed.txIRR.asV) * std.complex.expi(uniform(0, 1.0f, rnd) * 2 * PI);
        }

        /* RX IQ Mixer の設定 */
        {
            Random rnd = uniqueRandom(iTrial, "RXIQMixer");
            model.rxIQMixer.imbCoef = (1.0L / modelSeed.rxIRR.asV) * std.complex.expi(uniform(0, 1.0f, rnd) * 2 * PI);
        }

        /* チャネルの設定 */
        {
            Random rnd = uniqueRandom(iTrial, "Channel");
            model.channel.taps = 64;

            BoxMuller!Random gGen = BoxMuller!Random(rnd);
            Complex!real[] coefs;
            foreach(i; 0 .. model.channel.taps){
                // 64タップで40dB減衰
                auto db = -1 * (40.0 / 64.0) * i;
                coefs ~= cast(Complex!real)(gGen.front * 10.0L ^^ (db/20));
                gGen.popFront();
            }

            model.channel.impulseResponse = coefs;
        }

        /* キャンセラの設定 */
        {
            model.orthogonalizer.numOfTrainingSymbols = 10000;
            model.firFilter.taps = 64;

            if(methodName[0] == 'O')
                model.orthogonalizer.enabled = true;
            else
                model.orthogonalizer.enabled = false;

            // ベースバンド信号波形の出力
            model.outputWaveform = false;

            // BERの出力
            model.outputBER = true;

            if(methodName.canFind("DCM")) {
                model.learningSymbols = 10;
                model.learningCount = 3;
            } else {
                model.learningSymbols = modelSeed.numOfTrainingSymbols;
                model.learningCount = 1;
            }

            model.swappedSymbols = 100000;
            model.rlsAdapter.delta = modelSeed.rlsDelta;
            model.rlsAdapter.lambda = modelSeed.rlsLambda;
            model.nlmsAdapter.mu = modelSeed.nlmsMu;
            model.iterativeFreqSIC.iterations = modelSeed.iterNumOfIteration;
        }
    }

    return models;
}


string makeDirNameOfModelSeed(ModelSeed modelSeed)
{
    string dir;
    if(modelSeed.cancellerType.split("_")[0].endsWith("FHF"))
    {
        dir = "TXP%s_inr%s_snr%s_%s_SF%s_G%s_IRR%s_%s"
            .format(modelSeed.txPower,
                    modelSeed.INR,
                    modelSeed.SNR,
                    modelSeed.cancellerType,
                    modelSeed.lnaSmoothFactor,
                    modelSeed.bfsNoiseMargin,
                    modelSeed.txIRR,
                    modelSeed.numOfTrainingSymbols);
    }
    else
    {
        dir = "TXP%s_inr%s_snr%s_%s_SF%s_IRR%s_%s"
            .format(modelSeed.txPower,
                    modelSeed.INR,
                    modelSeed.SNR,
                    modelSeed.cancellerType,
                    modelSeed.lnaSmoothFactor,
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
    Model[] models = makeModels!methodName(nTrials, modelSeed);

    if(exists(buildPath(dir, "allResult.json"))) return;

    JSONValue[] resList;
    uint sumOfSuccFreq;
    JSONValue[] selectingRatioList;
    foreach(i, ref m; models) {
        auto res = mainImpl!methodName(m, i == 0 ? dir : null);

        if(saveAllRAWData) {
            res["model"] = (){
                static JSONValue cpx2JV(F)(Complex!F a) { return JSONValue(["re": a.re, "im": a.im]); }

                JSONValue jv = cast(JSONValue[string])null;
                jv["impulseResponse"] = m.channel.impulseResponse.map!cpx2JV.array();
                jv["txIQCoef"] = cpx2JV(m.txIQMixer.imbCoef);
                jv["rxIQCoef"] = cpx2JV(m.rxIQMixer.imbCoef);

                return jv;
            }();
        }
        resList ~= res;

        if(methodName.startsWith("SFHF") || methodName.startsWith("S1FHF") || methodName.startsWith("S2FHF")){
            auto cnt = res["filterSpec"]["selectingIsSuccess"].array.map!(a => a.type == JSON_TYPE.TRUE ? 1 : 0).sum();
            sumOfSuccFreq += cnt;
            selectingRatioList ~= JSONValue(cnt / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling));
        }
    }

    JSONValue jv = cast(JSONValue[string])null;
    jv["cancellations"] = resList.map!(a => a["cancellation_dB"]).array();
    jv["bers"] = resList.map!(a => a["ber"]).array();
    
    auto file = File(buildPath(dir, "allResult.json"), "w");
    file.write(jv.toPrettyString(JSONOptions.specialFloatLiterals));

    if(saveAllRAWData){
        file = File(buildPath(dir, "rawAllResult.json"), "w");
        file.write(JSONValue(resList).toPrettyString(JSONOptions.specialFloatLiterals));
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
