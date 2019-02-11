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
    // Voltage txPower = 23.dBm;
    // Voltage paIIP3 = 20.dBm;
    // Gain paGain = 27.dB;
    enum real txBackoff_dB = 10;
    enum real txPower_dBm = 23;
    enum real paGain_dB = 20;
    Voltage txPower = txPower_dBm.dBm;
    Voltage paIIP3 = (txPower_dBm - paGain_dB + 6 + txBackoff_dB).dBm;
    Gain paGain = paGain_dB.dB;

    /* LNA */
    uint lnaSmoothFactor = 1;

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
    uint iterNumOfIteration = 3;
    uint iterNumOfNewton = 2;

    /* measure only desired signal */
    bool onlyDesired = false;

    bool outputBER = false;
    bool outputEVM = false;
}



void mainJob()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    auto taskList = new MultiTaskList();

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


    enum numOfTrials = 11;

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
                                    // "PH_LS",
                                    // "PH_RLS",
                                    // "PH_RLS",
                                    // "OPH_LMS",
                                    
                                    // "WL_LS",
                                    // "L_LS",
                                    
                                    // "IterativeFreqSIC_X",
                                    // "SidelobeFwd_X",
                                    // "SidelobeInv_X",
                                    "SidelobeInv2_X",
            ))
    {
        bool[string] dirset;
        auto appender = uniqueTaskAppender(&mainForEachTrial!methodName);
        scope(exit){
            enforce(appender.length == dirset.length);
            taskList ~= appender;
        }

        /+
        /* change the number of iterations */
        static if(methodName == "SidelobeInv2_X")
        foreach(nIters; iota(1, 11))
        {
            ModelSeed modelSeed;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = 10;
            modelSeed.INR = 60.dB;
            modelSeed.outputBER = false;
            modelSeed.outputEVM = false;
            modelSeed.iterNumOfIteration = nIters;
            modelSeed.iterNumOfNewton = 10;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath("results_estimate_iters", dir ~ format("_%s", nIters));
            dirset[dir] = true;
            appender.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }


        /* change the number of newton's method loop */
        static if(methodName == "SidelobeInv2_X")
        foreach(newtonIters; iota(0, 11))
        {
            ModelSeed modelSeed;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = 10;
            modelSeed.INR = 60.dB;
            modelSeed.outputBER = false;
            modelSeed.outputEVM = false;
            modelSeed.iterNumOfIteration = 10;
            modelSeed.iterNumOfNewton = newtonIters;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath("results_newton_iters", dir ~ format("_%s", newtonIters));
            dirset[dir] = true;
            appender.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }


        // only desired signal
        static if(methodName == "PH_LS")
        foreach(snr; iota(0, 21, 1))
        {
            ModelSeed modelSeed;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = 10;
            modelSeed.INR = 20.dB;
            modelSeed.SNR = snr.dB;
            modelSeed.onlyDesired = true;
            modelSeed.outputBER = true;
            modelSeed.outputEVM = true;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath("results_ber", dir ~ "_onlyDesired");
            dirset[dir] = true;
            appender.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }


        // learning symbols vs (EVM / SIC / BER)
        foreach(inr; iota(50, 62, 2))
        foreach(learningSymbols; iota(2, 21))
        {
            ModelSeed modelSeed;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = learningSymbols;
            modelSeed.INR = inr.dB;
            modelSeed.SNR = 20.dB;
            modelSeed.outputBER = true;
            modelSeed.outputEVM = true;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath("results_ber", dir);
            dirset[dir] = true;
            appender.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }
        +/

        // INR vs (EVM / SIC / BER)
        foreach(inr; iota(20, 82, 5))
        {
            ModelSeed modelSeed;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = 20;
            modelSeed.INR = inr.dB;
            modelSeed.outputBER = false;
            modelSeed.outputEVM = false;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath("results_inr_vs_sic", dir);
            dirset[dir] = true;
            appender.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }

        /+
        // TXP vs (EVM. SIC /  BER)
        foreach(txp; iota(10, 32, 1)) {
            ModelSeed modelSeed;
            modelSeed.txPower = txp.dBm;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = 10;
            modelSeed.INR = ((txp - 23)+50).dB;
            modelSeed.txPower = txp.dBm;
            modelSeed.outputBER = true;
            modelSeed.outputEVM = true;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath("results_txp_vs_sic", dir);
            dirset[dir] = true;
            appender.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }
        +/
    }

    import std.stdio;

    // writefln("%s tasks will be submitted.", taskList.length);
    JobEnvironment env;
    tuthpc.taskqueue.run(taskList, env);
    // foreach(i; 0 .. taskList.length)
    //     taskList[i]();
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
            model.channel.taps = 8;

            BoxMuller!Random gGen = BoxMuller!Random(rnd);
            Complex!real[] coefs;
            foreach(i; 0 .. model.channel.taps){
                // tapsのタップ数で40dB減衰
                auto db = -1 * (40.0 / model.channel.taps) * i;
                coefs ~= cast(Complex!real)(gGen.front * 10.0L ^^ (db/20));
                gGen.popFront();
            }

            model.channel.impulseResponse = coefs;
        }

        /* キャンセラの設定 */
        {
            model.orthogonalizer.numOfTrainingSymbols = 10000;
            model.firFilter.taps = model.channel.taps;

            if(methodName[0] == 'O')
                model.orthogonalizer.enabled = true;
            else
                model.orthogonalizer.enabled = false;

            // ベースバンド信号波形の出力
            model.outputWaveform = false;

            // BERの出力
            model.outputBER = modelSeed.outputBER;
            model.outputEVM = modelSeed.outputEVM;

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
            model.iterativeFreqSIC.newtonIterations = modelSeed.iterNumOfNewton;
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
        {
            import core.memory;
            GC.collect();
        }
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
    if(modelSeed.outputBER) jv["bers"] = resList.map!(a => a["ber"]).array();
    if(modelSeed.outputEVM) jv["evms"] = resList.map!(a => a["evm"]).array();
    
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
