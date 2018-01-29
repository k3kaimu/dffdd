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

import dffdd.utils.unit;

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
    /* IQ Mixer */
    Gain txIRR = 25.dB;
    Gain rxIRR = 25.dB;

    /* PA */
    Voltage txPower = 23.dBm;
    Voltage paIIP3 = 21.8.dBm;
    Gain paGain = 28.5.dB;

    /* LNA */
    uint lnaSmoothFactor = 1;

    /* Other */
    Gain snr = 11.dB;
    Gain inr = 50.dB;
    uint numOfTrainingSymbols = 20;
    Gain gamma = 0.dB;
}



void mainJob()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    auto taskList = new MultiTaskList();

    auto setRLSLMSParam(string methodName, ref Model model)
    {
        if(methodName.startsWith("OPH")) {
            model.rlsAdapter.delta = 0.1;
            model.rlsAdapter.lambda = 1 - 1.5E-3;
            model.nlmsAdapter.mu = 0.33;
        } else {
            model.rlsAdapter.delta = 4E-7;
            model.rlsAdapter.lambda = 1;
            model.nlmsAdapter.mu = 1;
        }
    }


    enum numOfTrials = 100;

    // ADC&IQ&PA
    foreach(methodName; AliasSeq!(
                                    //"S2FHF_LS",
                                    //"S2FHF_RLS",
                                    //"S2FHF_LMS",
                                    ////
                                    //"FHF_LS",
                                    //"FHF_RLS",
                                    //"FHF_LMS",
                                    ////
                                    //"OPH_RLS",
                                    //"OPH_LS",
                                    //"OPH_LMS",
                                    //
                                    //"WL_LS",
                                    //
                                    //"L_LS",
                                    //
                                    "IterativeFreqSIC_X",
            ))
    {
        bool[string] dirset;
        auto appender = uniqueTaskAppender(&mainForEachTrial!methodName);
        scope(exit){
            enforce(appender.length == dirset.length);
            taskList ~= appender;
        }

        /+
        /// S2FHF: gamma vs cancellation
        static if(methodName == "S2FHF_LS")
        foreach(learningSymbols; iota(2, 22, 2))
        foreach(gamma; iota(-10, 12, 2))
        {
            auto md = makeModelAndDir!methodName(learningSymbols, 50, 23, 3, gamma, 25);
            auto dir = buildPath("results_S2FHFLS_gamma_vs_canc", md[1]);
            appender.append(md[0], dir, No.saveAllRAWData);
            dirset[dir] = true;
        }


        /// S2FHF: selection
        static if(methodName == "S2FHF_LS")
        foreach(inr; iota(20, 75, 5))
        foreach(sf; [1, 3])
        {
            auto md = makeModelAndDir!methodName(20, inr, 23, sf, -4, 25);
            auto dir = buildPath("results_S2FHFLS_selection", md[1]);
            appender.append(md[0], dir, Yes.saveAllRAWData);
            dirset[dir] = true;
        }


        /// Iterative: iter vs cancellation
        static if(methodName == "IterativeFreqSIC_X")
        foreach(learningSymbols; iota(2, 11, 1))
        foreach(inr; [10, 20, 30, 40, 50, 60, 70])
        foreach(iter; iota(1, 6))
        {
            auto md = makeModelAndDir!methodName(learningSymbols, inr, 23, 3, -4, 25);
            auto dir = buildPath("results_Iterative_iter_vs_canc", md[1]);
            dir ~= format("_%s", iter);
            md[0].iterativeFreqSIC.iterations = iter;
            appender.append(md[0], dir, No.saveAllRAWData);
            dirset[dir] = true;
        }


        /// inr vs cancellation
        static if(methodName.endsWith("_LS") || methodName == "IterativeFreqSIC_X")
        foreach(learningSymbols; [10, 20, 30, 100])
        foreach(inr; iota(20, 75, 5))
        foreach(sf; [1, 3])
        {
            auto md = makeModelAndDir!methodName(learningSymbols, inr, 23, sf, -4, 25);
            auto dir = buildPath("results_ALL_inr_vs_canc", md[1]);
            appender.append(md[0], dir, No.saveAllRAWData);
            dirset[dir] = true;
        }
        +/


        foreach(learningSymbols; iota(2, 21, 1))
        {
            auto md = makeModelAndDir!methodName(numOfTrials, learningSymbols, 50, 23, 3, -4, 25);
            auto dir = buildPath("results_ALL_iteration_vs_canc", md[1]);
            setRLSLMSParam(methodName, md[0]);
            appender.append(md[0], dir, No.saveAllRAWData);
            dirset[dir] = true;
        }


        ///// RLS parameter, lambda = 1, delta=???
        //static if(methodName.endsWith("_RLS"))
        //foreach(lambdaExp; [-5, -4, /*-2, -3, -4, -5, -6, -7, -8, -9*/])
        //foreach(lambdaDig; iota(10, 100)){
        //    auto md = makeModelAndDir!methodName(10, 50, 23, 3, -4, 25);
        //    auto dir = buildPath("results_RLSLMSParam", md[1]);
        //    dir ~= format("_%s_%s_deltaopt", lambdaDig, lambdaExp);
        //    with(md[0].rlsAdapter){
        //        lambda = 1 - (lambdaDig * (10.0L^^lambdaExp));
        //        if(methodName.startsWith("FHF"))
        //            delta = 4E-7;
        //        else
        //            delta = 0.1;
        //    }
        //    appender.append(md[0], dir, No.saveAllRAWData);
        //    dirset[dir] = true;
        //}

        //static if(methodName.endsWith("_LMS"))
        //foreach(muExp; [-2])
        //foreach(muDig; iota(21, 51, 1)){
        //    auto md = makeModelAndDir!methodName(10, 50, 23, 3, -4, 25);
        //    auto dir = buildPath("results_RLSLMSParam", md[1]);
        //    dir ~= format("_%s_%s", muDig, muExp);
        //    with(md[0].nlmsAdapter){
        //        mu = muDig * (10.0L^^muExp);
        //    }
        //    appender.append(md[0], dir, No.saveAllRAWData);
        //    dirset[dir] = true;
        //}
    }

    //writefln("%s tasks will be submitted.", taskList.length);
    JobEnvironment env;
    tuthpc.taskqueue.run(taskList, env);
}


Tuple!(Model[], string) makeModelAndDir(string methodName)(size_t numOfTrials, ModelSeed modelSeed)
{
    Model[] models;

    foreach(iTrial; 0 .. numOfTrials) {
        Model model;
        scope(success)
            models ~= model;

        model.SNR = modelSeed.snr;
        model.INR = modelSeed.inr;
        model.pa.TX_POWER = modelSeed.txPower;
        model.txIQMixer.IRR = modelSeed.txIRR;
        model.rxIQMixer.IRR = modelSeed.rxIRR;
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
        model.useSTXIQ2 = false;
        model.useSRXLN = true;
        model.useSRXIQ = true;
        model.useSRXQZ = true;


        /* PAの設定 */
        {
            model.pa.TX_POWER = modelSeed.txPower;
            model.pa.IIP3 = modelSeed.paIIP3.dBm;
            model.pa.GAIN = modelSeed.paGain.dB;
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
            model.channel.isCoaxialCable = false;
        }

        /* キャンセラの設定 */
        {
            model.orthogonalizer.numOfTrainingSymbols = 10000;

            if(methodName[0] == 'O')
                model.orthogonalizer.enabled = true;
            else
                model.orthogonalizer.enabled = false;

            // ベースバンド信号波形の出力
            model.outputWaveform = false;

            model.channel.taps = 64;
            model.firFilter.taps = 64;

            model.outputBER = false;

            if(methodName.canFind("DCM")) {
                model.learningSymbols = 10;
                model.learningCount = 3;
            } else {
                model.learningSymbols = learningSymbols;
                model.learningCount = 1;
            }

            if(methodName.split("_")[0].endsWith("FHF") || methodName == "IterativeFreqSIC_X") {
                model.swappedSymbols = 100000;
                model.rlsAdapter.delta = 4E-7;
                model.rlsAdapter.lambda = 1;
                model.nlmsAdapter.mu = 1;
            } else {
                model.swappedSymbols = 0;
                model.rlsAdapter.delta = 0.1;
                model.rlsAdapter.lambda = 1;
                model.nlmsAdapter.mu = 0.33;
            }
        }
    }

  static if(methodName.split("_")[0].endsWith("FHF"))
    string dir = "TXP%s_inr%s_%s_SF%s_G%s_IRR%s_%s".format(model.pa.TX_POWER, model.INR, methodName, model.lna.smoothFactor, model.basisFuncsSelection.noiseMargin, irr, learningSymbols);
  else
    string dir = "TXP%s_inr%s_%s_SF%s_IRR%s_%s".format(model.pa.TX_POWER, model.INR, methodName, model.lna.smoothFactor, irr, learningSymbols);

    return typeof(return)(model, dir);
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


void mainForEachTrial(string methodName)(Model m, string dir, Flag!"saveAllRAWData" saveAllRAWData = No.saveAllRAWData)
{
    if(exists(buildPath(dir, "allResult.json"))) return;

    JSONValue[] resList;

    // 最初の一回は普通にやる
    resList ~= mainImpl!methodName(m, dir);

    // writeln(mainImpl!methodName(m, dir)["training_symbols_per_second"]);
    enum K = 100;    // 試行回数
    uint sumOfSuccFreq;
    JSONValue[] selectingRatioList;
    foreach(j; 0 .. K){
        m.rndSeed += 100;   // seed値を100ずつ足していく
        auto res = mainImpl!methodName(m, null);
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
            JSONValue[string] vv;
            vv["selectingRatio"] = sumOfSuccFreq / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * K);
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
                foreach(f; 0 .. m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling){
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

                foreach(f; 0 .. m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling)
                    foreach(b; used[f].array)
                        if(b.type == JSON_TYPE.TRUE)
                            ++selectedBF;
            }
            vv["selectedBFRatio"] = selectedBF / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length * distDim);
            vv["avgSelectedBF"] = selectedBF / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);
            vv["failedRatio"] = failedCNT / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);
            vv["matchRatio"] = matchCNT / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);
            vv["overRatio"] = overCNT / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);

            vv["failedPLQRatio"] = failedCNTPLQ / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);
            vv["matchPLQRatio"] = matchCNTPLQ / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);
            vv["overPLQRatio"] = overCNTPLQ / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * resList.length);

            immutable size_t nfft = m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling;
            immutable size_t nsym = (m.ofdm.numOfFFT + m.ofdm.numOfCP) * m.ofdm.scaleOfUpSampling;
            vv["idealCostRLS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + (idealCostRLS / resList.length))/nsym;
            vv["idealCostLMS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + (idealCostLMS / resList.length))/nsym;
            vv["actualCostRLS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + (actualCostRLS / resList.length))/nsym;
            vv["actualCostLMS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + (actualCostLMS / resList.length))/nsym;
            vv["nonSelCostRLS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + 4*(distDim^^2 + distDim) * nfft)/nsym;
            vv["nonSelCostLMS"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + 2 * distDim * nfft)/nsym;
            vv["nonSelCostRLS2"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + 4*(2^^2 + 2) * nfft)/nsym;
            vv["nonSelCostLMS2"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + 2 * 2 * nfft)/nsym;
            vv["nonSelCostRLS3"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + 4*(3^^2 + 2) * nfft)/nsym;
            vv["nonSelCostLMS3"] = (0.25 * (distDim + 2) * nfft * log2(nfft) + 3 * 2 * nfft)/nsym;

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
