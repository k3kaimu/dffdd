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


    enum numOfTrials = 5;

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
                                    
                                    "IterativeFreqSIC_X",
                                    "Sidelobe_X",
            ))
    {
        bool[string] dirset;
        auto appender = uniqueTaskAppender(&mainForEachTrial!methodName);
        scope(exit){
            enforce(appender.length == dirset.length);
            taskList ~= appender;
        }

        
        // /// S2FHF: gamma vs cancellation
        // static if(methodName == "S2FHF_LS")
        // foreach(learningSymbols; iota(2, 22, 2))
        // foreach(gamma; iota(-10, 12, 2))
        // {
        //     ModelSeed modelSeed;
        //     modelSeed.cancellerType = methodName;
        //     modelSeed.numOfTrainingSymbols = learningSymbols;
        //     modelSeed.bfsNoiseMargin = gamma.dB;

        //     auto dir = makeDirNameOfModelSeed(modelSeed);
        //     dir = buildPath("results_S2FHFLS_gamma_vs_canc", dir);
        //     dirset[dir] = true;

        //     appender.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        // }


        // /// S2FHF: selection
        // static if(methodName == "S2FHF_LS")
        // foreach(inr; iota(20, 75, 5))
        // foreach(sf; [1, 3])
        // {
        //     ModelSeed modelSeed;
        //     modelSeed.cancellerType = methodName;
        //     modelSeed.numOfTrainingSymbols = 20;
        //     modelSeed.INR = inr.dB;
        //     modelSeed.lnaSmoothFactor = sf;

        //     auto dir = makeDirNameOfModelSeed(modelSeed);
        //     dir = buildPath("results_S2FHFLS_selection", dir);
        //     dirset[dir] = true;

        //     appender.append(numOfTrials, modelSeed, dir, Yes.saveAllRAWData);
        // }


        // /// Iterative: iter vs cancellation
        // static if(methodName == "IterativeFreqSIC_X")
        // foreach(learningSymbols; iota(2, 11, 1))
        // foreach(inr; [10, 20, 30, 40, 50, 60, 70])
        // foreach(iter; iota(1, 5))
        // {
        //     ModelSeed modelSeed;
        //     modelSeed.cancellerType = methodName;
        //     modelSeed.numOfTrainingSymbols = learningSymbols;
        //     modelSeed.INR = inr.dB;
        //     modelSeed.iterNumOfIteration = iter;

        //     auto dir = makeDirNameOfModelSeed(modelSeed);
        //     dir ~= format("_%s", iter);
        //     dir = buildPath("results_Iterative_iter_vs_canc", dir);
        //     dirset[dir] = true;

        //     appender.append(numOfTrials, modelSeed, dir, Yes.saveAllRAWData);
        // }


        /// inr vs cancellation
        // static if(methodName.endsWith("_LS") || methodName == "IterativeFreqSIC_X")
        foreach(learningSymbols; [50])
        foreach(inr; iota(20, 105, 1))
        foreach(sf; [1, 3])
        {
            ModelSeed modelSeed;
            modelSeed.cancellerType = methodName;
            modelSeed.numOfTrainingSymbols = learningSymbols;
            modelSeed.INR = inr.dB;
            modelSeed.lnaSmoothFactor = sf;

            auto dir = makeDirNameOfModelSeed(modelSeed);
            dir = buildPath("results_ALL_inr_vs_canc", dir);
            dirset[dir] = true;

            appender.append(numOfTrials, modelSeed, dir, No.saveAllRAWData);
        }


        // foreach(learningSymbols; iota(2, 21, 1))
        // {
        //     ModelSeed modelSeed;
        //     modelSeed.cancellerType = methodName;
        //     modelSeed.numOfTrainingSymbols = learningSymbols;
        //     setRLSLMSParam(modelSeed);
        //     auto dirName = makeDirNameOfModelSeed(modelSeed);
        //     dirName = buildPath("results_ALL_iteration_vs_canc", dirName);

        //     appender.append(numOfTrials, modelSeed, dirName, No.saveAllRAWData);
        //     dirset[dirName] = true;
        // }


        // // RLS parameter, lambda = 1, delta=???
        // static if(methodName.endsWith("_RLS"))
        // foreach(lambdaExp; [-5, -4, /*-2, -3, -4, -5, -6, -7, -8, -9*/])
        // foreach(lambdaDig; iota(10, 100)) {
        //     ModelSeed modelSeed;
        //     modelSeed.cancellerType = methodName;
        //     modelSeed.numOfTrainingSymbols = 10;
        //     modelSeed.rlsLambda = 1 - (lambdaDig * (10.0L^^lambdaExp));
            
        //     if(modelSeed.cancellerType.startsWith("FHF")) {
        //         modelSeed.rlsDelta = 4E-7;
        //     }else{
        //         modelSeed.rlsDelta = 0.1;
        //     }

        //     auto dirName = makeDirNameOfModelSeed(modelSeed);
        //     dirName ~= format("_%s_%s_deltaopt", lambdaDig, lambdaExp);
        //     dirName = buildPath("results_RLSLMSParam", dirName);
        //     dirset[dirName] = true;

        //     appender.append(numOfTrials, modelSeed, dirName, No.saveAllRAWData);
        // }


        // static if(methodName.endsWith("_LMS"))
        // foreach(muExp; [-2])
        // foreach(muDig; iota(21, 51, 1)){
        //     ModelSeed modelSeed;
        //     modelSeed.cancellerType = methodName;
        //     modelSeed.numOfTrainingSymbols = 10;
        //     modelSeed.nlmsMu = muDig * (10.0L^^muExp);

        //     string dirName = makeDirNameOfModelSeed(modelSeed);
        //     dirName ~= format("_%s_%s", muDig, muExp);
        //     dirName = buildPath("results_RLSLMSParam", dirName);

        //     appender.append(numOfTrials, modelSeed, dirName, No.saveAllRAWData);
        //     dirset[dirName] = true;
        // }
    }

    //writefln("%s tasks will be submitted.", taskList.length);
    JobEnvironment env;
    tuthpc.taskqueue.run(taskList, env);
}


Model[] makeModels(string methodName)(size_t numOfTrials, ModelSeed modelSeed)
{
    Model[] models;

    foreach(iTrial; 0 .. numOfTrials) {
        Model model;
        scope(success)
            models ~= model;

        model.rndSeed = iTrial + 0xDEAD + 0xBEAF;
        model.SNR = modelSeed.SNR;
        model.INR = modelSeed.INR;
        // model.pa.TX_POWER = modelSeed.txPower;
        // model.txIQMixer.IRR = modelSeed.txIRR;
        // model.rxIQMixer.IRR = modelSeed.rxIRR;
        model.basisFuncsSelection.noiseMargin = modelSeed.gamma;
        model.lna.smoothFactor = modelSeed.lnaSmoothFactor;

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
            model.txIQMixer.imbCoef = (1.0L / modelSeed.txIRR.gain) * std.complex.expi(uniform(0, 1.0f, rnd) * 2 * PI);
        }

        /* RX IQ Mixer の設定 */
        {
            Random rnd = uniqueRandom(iTrial, "RXIQMixer");
            model.rxIQMixer.imbCoef = (1.0L / modelSeed.rxIRR.gain) * std.complex.expi(uniform(0, 1.0f, rnd) * 2 * PI);
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
            model.outputBER = false;

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
        dir = "TXP%s_inr%s_%s_SF%s_G%s_IRR%s_%s"
            .format(modelSeed.txPower,
                    modelSeed.INR,
                    modelSeed.cancellerType,
                    modelSeed.lnaSmoothFactor,
                    modelSeed.bfsNoiseMargin,
                    modelSeed.txIRR,
                    modelSeed.numOfTrainingSymbols);
    }
    else
    {
        dir = "TXP%s_inr%s_%s_SF%s_IRR%s_%s"
            .format(modelSeed.txPower,
                    modelSeed.INR,
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

    if(methodName.startsWith("SFHF") || methodName.startsWith("S1FHF") || methodName.startsWith("S2FHF"))
        jv["selecting"] = (){
            immutable nFFT = models[0].ofdm.numOfFFT * models[0].ofdm.scaleOfUpSampling;

            JSONValue[string] vv;
            vv["selectingRatio"] = sumOfSuccFreq / cast(float)(nFFT * nTrials);
            vv["selectingRatioList"] = selectingRatioList;

            size_t selectedBF;
            size_t failedCNT, failedCNTPLQ;
            size_t matchCNT, matchCNTPLQ;
            size_t overCNT, overCNTPLQ;
            real idealCostRLS = 0;
            real idealCostLMS = 0;
            real actualCostRLS = 0;
            real actualCostLMS = 0;
            size_t distDim = 0;
            foreach(k, ref ev; resList){
                auto reqs = ev["filterSpec"]["actualRequiredBasisFuncs"].array;
                auto used = ev["filterSpec"]["selectedBasisFuncs"].array;
                foreach(f; 0 .. nFFT){
                    auto r = reqs[f].array;
                    auto u = used[f].array;

                    if(distDim == 0) distDim = r.length;

                    // idealCostに追加する
                    size_t reqsP, usedP;
                    foreach(p; 0 .. r.length){
                        if(r[p].type == JSON_TYPE.TRUE) ++reqsP;
                        if(u[p].type == JSON_TYPE.TRUE) ++usedP;
                    }
                    idealCostRLS += 4*(reqsP^^2) + 4*reqsP;
                    idealCostLMS += 2*reqsP;
                    actualCostRLS += 4*(usedP^^2) + 4*usedP;
                    actualCostLMS += 2*usedP;

                    // 推定ミスした周波数のカウント
                    foreach(p; 0 .. r.length) if(r[p].type == JSON_TYPE.TRUE && u[p].type == JSON_TYPE.FALSE) { ++failedCNT; break; }

                    // 過大評価した周波数のカウント
                    foreach(p; 0 .. r.length) if(r[p].type == JSON_TYPE.FALSE && u[p].type == JSON_TYPE.TRUE){ ++overCNT; break; }
                    
                    // 完全推定できた周波数のカウント
                    foreach(p; 0 .. r.length){
                        if(r[p].type != u[p].type) break;
                        if(p == r.length - 1) ++matchCNT;
                    }

                    // p>qについて
                    // 推定ミスした周波数のカウント
                    foreach(p; 0 .. r.length) if(CompleteDistorter!().indexOfConjugated(p) >= p) if(r[p].type == JSON_TYPE.TRUE && u[p].type == JSON_TYPE.FALSE) { ++failedCNTPLQ; break; }

                    // 過大評価した周波数のカウント
                    foreach(p; 0 .. r.length) if(CompleteDistorter!().indexOfConjugated(p) >= p) if(r[p].type == JSON_TYPE.FALSE && u[p].type == JSON_TYPE.TRUE){ ++overCNTPLQ; break; }
                    
                    // 完全推定できた周波数のカウント
                    bool bBreaked = false;
                    foreach(p; 0 .. r.length) if(CompleteDistorter!().indexOfConjugated(p) >= p) {
                        if(r[p].type != u[p].type){
                            bBreaked = true;
                            break;
                        }
                    }
                    if(!bBreaked) ++matchCNTPLQ;
                }

                foreach(f; 0 .. nFFT)
                    foreach(b; used[f].array)
                        if(b.type == JSON_TYPE.TRUE)
                            ++selectedBF;
            }
            vv["selectedBFRatio"] = selectedBF / cast(float)(nFFT * nTrials * distDim);
            vv["avgSelectedBF"] = selectedBF / cast(float)(nFFT * nTrials);
            vv["failedRatio"] = failedCNT / cast(float)(nFFT * nTrials);
            vv["matchRatio"] = matchCNT / cast(float)(nFFT * nTrials);
            vv["overRatio"] = overCNT / cast(float)(nFFT * nTrials);

            vv["failedPLQRatio"] = failedCNTPLQ / cast(float)(nFFT * nTrials);
            vv["matchPLQRatio"] = matchCNTPLQ / cast(float)(nFFT * nTrials);
            vv["overPLQRatio"] = overCNTPLQ / cast(float)(nFFT * nTrials);

            immutable size_t nsym = models[0].ofdm.numOfSamplesOf1Symbol;
            vv["idealCostRLS"] = (0.25 * (distDim + 2) * nFFT * log2(nFFT) + (idealCostRLS / nTrials))/nsym;
            vv["idealCostLMS"] = (0.25 * (distDim + 2) * nFFT * log2(nFFT) + (idealCostLMS / nTrials))/nsym;
            vv["actualCostRLS"] = (0.25 * (distDim + 2) * nFFT * log2(nFFT) + (actualCostRLS / nTrials))/nsym;
            vv["actualCostLMS"] = (0.25 * (distDim + 2) * nFFT * log2(nFFT) + (actualCostLMS / nTrials))/nsym;
            vv["nonSelCostRLS"] = (0.25 * (distDim + 2) * nFFT * log2(nFFT) + 4*(distDim^^2 + distDim) * nFFT)/nsym;
            vv["nonSelCostLMS"] = (0.25 * (distDim + 2) * nFFT * log2(nFFT) + 2 * distDim * nFFT)/nsym;
            vv["nonSelCostRLS2"] = (0.25 * (distDim + 2) * nFFT * log2(nFFT) + 4*(2^^2 + 2) * nFFT)/nsym;
            vv["nonSelCostLMS2"] = (0.25 * (distDim + 2) * nFFT * log2(nFFT) + 2 * 2 * nFFT)/nsym;
            vv["nonSelCostRLS3"] = (0.25 * (distDim + 2) * nFFT * log2(nFFT) + 4*(3^^2 + 2) * nFFT)/nsym;
            vv["nonSelCostLMS3"] = (0.25 * (distDim + 2) * nFFT * log2(nFFT) + 3 * 2 * nFFT)/nsym;

            return vv;
        }();
    
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
