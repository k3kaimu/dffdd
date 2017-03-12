module app;


import core.thread;

import std.random;
import std.algorithm;
import std.range;
import std.conv;
import std.complex;
import std.stdio;
//import msgpackrpc;
import std.datetime;
import std.math;
import std.functional;
import std.mathspecial;
import std.format;
import std.parallelism;
import std.path;
import std.file : mkdirRecurse;
import std.meta;
import std.numeric;
import std.exception;
import std.json;
import std.variant;


import carbon.math : nextPowOf2;

import dffdd.blockdiagram.utils;
import dffdd.utils.fft;
import dffdd.blockdiagram.amplifier;
import dffdd.utils.unit;
import dffdd.blockdiagram.noise;
import dffdd.blockdiagram.decimator;
import dffdd.blockdiagram.filter;
import dffdd.blockdiagram.adder;
import dffdd.blockdiagram.mod.ofdm;
import dffdd.blockdiagram.iqmixer;
import dffdd.blockdiagram.quantizer;
//import dffdd.utils.msgpackrpc;


import std.stdio;
import std.algorithm;
import std.typecons;
import std.range;
import std.random;

//import carbon.range;
import dranges.range;
import dranges.algorithm;

import dffdd.mod.primitives;
import dffdd.mod.qam;
import dffdd.mod.ofdm;

import dffdd.dsp.statistics;


import constant;
import models;
import snippet;
import sigsim;


real qfunc(real x) { return 0.5 * erfc(x / SQRT2); }


auto psdSaveTo(R)(R r, string filename, string resultDir, size_t dropSize, Model model)
if(isForwardRange!R)
{
    alias E = ElementType!R;
    return r.loggingTo(makeSpectrumAnalyzer!E(filename, resultDir, dropSize, model)).toWrappedRange;
}


auto makeSpectrumAnalyzer(C)(string filename, string resultDir, size_t dropSize, Model model)
{
    if(resultDir !is null) {
        return makeInstrument!C(delegate void(FiberRange!C r){
            r.drop(dropSize).writePSD(File(buildPath(resultDir, filename), "w"), model.samplingFreq, 1024);
        });
    }else{
        return makeInstrument!C(delegate void(FiberRange!C r){});
    }
}


auto makeCancellationProbe(C)(float* dst, string filename, string resultDir, size_t dropSize, Model model)
{
    alias Tup = Tuple!(C, C);

    return makeInstrument!Tup(delegate void(FiberRange!Tup r){
        auto rd = r.drop(dropSize);
        auto ratio = rd.calculateSIC(model.samplingFreq, 1024, model.ofdm.numOfFFT, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling);
        if(dst !is null) *dst = 10*log10(ratio);

        if(resultDir !is null) {
            auto file = File(buildPath(resultDir, filename), "w");
            file.writefln("%s", 10*log10(ratio));
        }
    });
}


auto makeCancellationIterationProbe(C)(string filename, string resultDir, size_t dropSize, Model model)
{
    alias Tup = Tuple!(C, C);

    if(resultDir !is null){
        return makeInstrument!Tup(delegate void(FiberRange!Tup r){
            auto rd = r.drop(dropSize);
            auto file = File(buildPath(resultDir, filename), "w");

            size_t cnt;
            while(!rd.empty){
                auto ratio = rd.calculateSIC(model.samplingFreq, 1024, model.ofdm.numOfFFT, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling, 1);
                //if(ratio.isNaN && rd.empty) return;
                if(ratio.isNaN) ratio = 0;

                file.writefln("%s,%s", cnt*1024, 10*log10(ratio));
                file.flush();
                ++cnt;
            }
        });
    }else{
        return makeInstrument!Tup(delegate void(FiberRange!Tup r) {});
    }
}


auto makeWaveformProbe(C)(string filename, string resultDir, size_t dropSize, Model model)
{
    if(resultDir !is null) {
        return makeInstrument!C(delegate void(FiberRange!C r){
            File datfile = File(buildPath(resultDir, filename.stripExtension ~ ".dat"), "w");
            File csvfile = File(buildPath(resultDir, filename.stripExtension ~ ".csv"), "w");
            r = r.drop(dropSize);

            float[2][] chunk = new float[2][](1024);
            while(!r.empty){
                chunk.length = 0;

                foreach(i; 0 .. 1024){
                    if(r.empty) break;
                    auto frnt = r.front;
                    chunk ~= [frnt.re, frnt.im];
                    r.popFront();
                }

                datfile.rawWrite(chunk);
                foreach(e; chunk)
                    csvfile.writefln("%s,%s,", e[0], e[1]);
            }
        });
    }else{
        return makeInstrument!C(delegate void(FiberRange!C r) {});
    }
}


JSONValue mainImpl(string filterType)(Model model, string resultDir = null)
{
    JSONValue infoResult = ["type": filterType];

    immutable bool* alwaysFalsePointer = new bool(false);
    immutable bool* alwaysTruePointer = new bool(true);

    immutable ofdmModSignalPower = (){
        auto _modOFDMTest = modOFDM(model);
        return randomBits(1, model).connectToModulator(_modOFDMTest, alwaysFalsePointer, model).measurePower(1024*1024);
    }();

    if(resultDir !is null)
        mkdirRecurse(resultDir);

    auto signals = makeSimulatedSignals(model, resultDir);
    signals.trainAGC();

    // フィルタの学習
  enum string filterStructure = filterType.split("_")[0];
  enum string filterOptimizer = filterType.split("_")[1];

  enum bool isOrthogonalized = filterStructure[0] == 'O';

  static if(filterStructure.endsWith("PHDCM"))
    auto filter = makeParallelHammersteinWithDCMethodFilter!isOrthogonalized(modOFDM(model), model);
  else static if(filterStructure.endsWith("ARPH"))
    auto filter = makeAliasRemovableParallelHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("PH"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("CH"))
    auto filter = makeCascadeHammersteinFilter!(filterOptimizer)(modOFDM(model), model);
//   else static if(filterStructure.endsWith("CWLH"))
//     auto filter = makeCascadeWLHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
//   else static if(filterStructure.endsWith("CWL1H"))
//     auto filter = makeCascadeWL1HammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("DCMFHF"))
  {
      static assert(!isOrthogonalized);
      enum string filterOption = filterStructure[0 .. $-6];

      static assert(filterOption.canFind('1') || filterOption.canFind('2'));
      static assert(filterOption.canFind('P') || filterOption.canFind('C'));

      enum size_t type = filterOption.canFind('1') ? 1 : 2;
      enum Flag!"isParallel" isParallel = filterOption.canFind('P') ? Yes.isParallel : No.isParallel;

      auto filter = makeFrequencyDCMHammersteinFilter!(type, isParallel, filterOptimizer)(model);
  }
  else static if(filterStructure.endsWith("CFHF"))
    static assert(0); // auto filter = makeFrequencyCascadeHammersteinFilter!(true, filterOptimizer)(model);
  else static if(filterStructure.endsWith("FHF")){
    static assert(!isOrthogonalized);
    auto filter = makeFrequencyHammersteinFilter!(filterOptimizer)(model, filterStructure.endsWith("SFHF"));
  }else static if(filterStructure.endsWith("WL"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, 1)(modOFDM(model), model);
  else static if(filterStructure.endsWith("L"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, 1, true)(modOFDM(model), model);
  else
    static assert("Cannot identify filter model.");

    // フィルタの学習
    auto recvs = new Complex!float[model.blockSize],
         refrs = new Complex!float[model.blockSize],
         outps = new Complex!float[model.blockSize];

    {
        assert(model.rndSeed != 893);

        filter.preLearning(model, (Model m){
            // m.rndSeed = 893;

            auto ss = makeSimulatedSignals(m);
            ss.trainAGC();

            return ss;
        });
    }

    {
        //received.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);
        //txReplica.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);

        auto fftObj = new Fft(model.blockSize.nextPowOf2(0));

        //File powerFile = File(buildPath(resultDir, "errorout_long.csv"), "w");
        auto sicIterValue = makeCancellationIterationProbe!(Complex!float)("errorout_long.csv", resultDir, 0, model);
        auto waveTransmit = makeWaveformProbe!(Complex!float)("waveform_transmit_init.csv", resultDir, 0, model);
        auto waveBeforeSIC = makeWaveformProbe!(Complex!float)("waveform_beforeSIC_init.csv", resultDir, 0, model);
        auto waveAfterSIC = makeWaveformProbe!(Complex!float)("waveform_afterSIC_init.csv", resultDir, 0, model);
        auto waveFilterOutput = makeWaveformProbe!(Complex!float)("waveform_filterOutput_init.csv", resultDir, 0, model);

        StopWatch sw;
        foreach(blockIdx; 0 .. model.numOfFilterTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol / model.blockSize)
        {
            static if(filterStructure.endsWith("FHF"))
            {
                if(blockIdx >= model.swappedSymbols * model.ofdm.numOfSamplesOf1Symbol / model.blockSize)
                    signals.fillBuffer!(["txBaseband", "receivedSI"])(refrs, recvs);
                else
                    signals.fillBuffer!(["txBasebandSWP", "receivedSISWP"])(refrs, recvs);
            }
            else
                signals.fillBuffer!(["txBaseband", "receivedSI"])(refrs, recvs);

            sw.start();
            filter.apply!(Yes.learning)(refrs, recvs, outps);
            sw.stop();
            //filter.apply!false(refrs, recvs, outps);

            if(model.outputWaveform){
                foreach(i; 0 .. model.blockSize){
                    .put(waveTransmit, refrs[i]);
                    .put(waveBeforeSIC, recvs[i]);
                    .put(waveAfterSIC, outps[i]);
                    .put(waveFilterOutput, recvs[i] - outps[i]);
                }
            }

            if(blockIdx == 0 && resultDir !is null)
            {
                File outFile = File(buildPath(resultDir, "errorout_start.csv"), "w");
                foreach(i; 0 .. model.blockSize){
                    if(i > 20000) break;
                    auto pr = recvs[i].re^^2 + recvs[i].im^^2,
                         po = outps[i].re^^2 + outps[i].im^^2,
                         c = -10*log10(po / pr);

                    outFile.writefln("%s,%s,%s,%s,", i, pr, po, c);
                }
            }

            .put(sicIterValue, recvs.zip(outps));
        }

        infoResult["training_symbols_per_second"] = model.numOfFilterTrainingSymbols / ((cast(real)sw.peek.usecs) / 1_000_000);
    }


    {
        float sicv;

        auto psdBeforePSD = makeSpectrumAnalyzer!(Complex!float)("psd_beforeSIC.csv", resultDir, 0, model);
        auto psdAfterSIC = makeSpectrumAnalyzer!(Complex!float)("psd_afterSIC.csv", resultDir, 0, model);
        auto foutPSD = makeSpectrumAnalyzer!(Complex!float)("psd_filter_output.csv", resultDir, 0, model);
        auto sicValue = makeCancellationProbe!(Complex!float)(&sicv, "cancellation_value.csv", resultDir, 0, model);

        StopWatch sw;
        foreach(blockIdxo; 0 .. 1024)
        {
            signals.fillBuffer!(["txBaseband", "receivedSI"])(refrs, recvs);

            sw.start();
            filter.apply!(No.learning)(refrs, recvs, outps);
            sw.stop();

            // blockIdxoが10以上になるまで，慣らし運転する
            if(blockIdxo > 10){
                foreach(i; 0 .. model.blockSize)
                    refrs[i] = recvs[i] - outps[i];

                .put(psdBeforePSD, recvs);
                .put(psdAfterSIC, outps);
                .put(foutPSD, refrs);
                .put(sicValue, recvs.zip(outps));
            }
        }

        infoResult["canceling_symbols_per_second"] = 1024 * model.blockSize / model.ofdm.numOfSamplesOf1Symbol / ((cast(real)sw.peek.usecs) / 1_000_000);
        infoResult["cancellation_dB"] = sicv;
    }

    //received.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);
    //txReplica.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);


    if(model.outputBER)
    {
        size_t inpSize, outSize;
        {
            auto mod = modOFDM(model);
            inpSize = mod.symInputLength;
            outSize = mod.symOutputLength;
        }

        // 復調器の学習に必要
        immutable numOfTrainingBits = model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol / outSize * inpSize;

        auto refBits = signals.desiredBaseband.save.connectToDemodulator(modOFDM(model), model);

        real berResult = -1;
        auto berCounter = makeInstrument(delegate void (FiberRange!(Complex!float) r){
            auto bits = r
            .connectTo!PowerControlAmplifier(ofdmModSignalPower.sqrt.V)
            .connectToDemodulator(modOFDM(model), model);

            // 慣らし運転
            bits.popFrontN(numOfTrainingBits);
            refBits.popFrontN(numOfTrainingBits);

            berResult = measureBER(bits, refBits, model.berCounter.totalBits);
        });

        auto rcvAfterSICPSD = makeSpectrumAnalyzer!(Complex!float)("psd_rcv_afterSIC.csv", resultDir, 0, model);
        auto rcvBeforeSICPSD = makeSpectrumAnalyzer!(Complex!float)("psd_rcv_beforeIC.csv", resultDir, 0, model);

        size_t loopCount;
        while(berResult == -1){
            ++loopCount;

            signals.fillBuffer!(["txBaseband", "received"])(refrs, recvs);

            if(model.withSIC && model.withSI)
                filter.apply!(No.learning)(refrs, recvs, outps);
            else
                outps[] = recvs[];

            foreach(e; outps){
                berCounter.put(Complex!float(e.re, e.im));

                // loopCountが10以上になるまで，慣らし運転する
                if(loopCount > 10)
                    rcvAfterSICPSD.put(e);
            }

            // loopCountが10以上になるまで，慣らし運転する
            if(loopCount > 10)
                foreach(e; recvs)
                    rcvBeforeSICPSD.put(e);
        }

        if(resultDir !is null){
            File file = File(buildPath(resultDir, "ber.csv"), "w");
            file.writeln(berResult);
        }

        infoResult["ber"] = berResult;

        // writefln("%s,%s,%2.0f,%2.0f,%s", model.withSI ? 1 : 0, model.withSIC ? 1 : 0, model.SNR.dB, model.INR.dB, berResult);
    }

    // ノイズ電力
    auto noisePSD = makeSpectrumAnalyzer!(Complex!float)("psd_noise_floor.csv", resultDir, 0, model);
    .put(noisePSD, signals.noise.save.map!(a => Complex!float(a)).take(64*1024));


    static if(is(typeof((){ filter.saveInfoToDir(""); })))
    {
        if(resultDir !is null)
            filter.saveInfoToDir(resultDir);
    }

    static if(is(typeof((){ infoResult["filterSpec"] = filter.info; })))
    {
        infoResult["filterSpec"] = filter.info;
    }

    if(resultDir !is null)
        std.file.write(buildPath(resultDir, "info.json"), infoResult.toPrettyString());

    return infoResult;
}


void mainJob()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    auto taskList = new MultiTaskList();

    // ADC&IQ&PA
    foreach(methodName; AliasSeq!(/*"FHF_LS",*/ /*"OPH_LS",*/ "SFHF_RLS"/*, "C1DCMFHF_LS", "C2DCMFHF_LS", "P1DCMFHF_LS", "P2DCMFHF_LS",*/
                // "FHF_LMS", "FHF_LS", "OPH_LS", "OPH_RLS", "OPH_LMS", "OCH_LS", "OCH_RLS", "OCH_LMS", "WL_LS", "WL_RLS", "WL_LMS", "L_LS", "L_RLS", "L_LMS" /*"FHF", "PH"*//*, "OPH", "OPHDCM", "OCH", "WL", "L",*/ /*"OPHDCM"*/
            ))
        foreach(learningSymbols; iota(60, 65, 5)) foreach(orthTrainingSymbols; [10000])
        {
            Model[] models;
            string[] dirs;

            foreach(inr; iota(20, 55, 5)) foreach(txp; iota(15, 18, 3))
            foreach(gamma; iota(0, 8, 2))
            foreach(beta; iota(0, 30, 5))
            {
                Model model;
                model.SNR = 11.dB;
                model.INR = inr.dB;
                model.pa.TX_POWER = txp.dBm;

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

                model.orthogonalizer.numOfTrainingSymbols = orthTrainingSymbols;

                if(methodName[0] == 'O')
                    model.orthogonalizer.enabled = true;
                else
                    model.orthogonalizer.enabled = false;

                // ベースバンド信号波形の出力
                model.outputWaveform = false;

                model.channel.taps = 64;
                model.firFilter.taps = 64;

                model.outputBER = false;

              static if(methodName.canFind("DCM"))
              {
                model.learningSymbols = 10;
                model.learningCount = 3;
              }
              else
              {
                model.learningSymbols = learningSymbols;
                model.learningCount = 1;
              }

              static if(methodName.split("_")[0].endsWith("FHF"))
              {
                model.swappedSymbols = 100000;
                model.numOfFilterTrainingSymbols = 1000;
              }
              else
              {
                model.swappedSymbols = 0;
                model.numOfFilterTrainingSymbols = 1000;
              }

                model.basisFuncsSelection.imageMargin = (-beta).dB;
                model.basisFuncsSelection.noiseMargin = gamma.dB;

                models ~= model;

              static if(methodName.endsWith("LMS") || methodName.endsWith("RLS"))
                model.numOfFilterTrainingSymbols = 100;

              static if(methodName.split("_")[0].endsWith("FHF"))
                dirs ~= "TXP%s_inr%s_%s_B%s_G%s".format(model.pa.TX_POWER, model.INR, methodName, model.basisFuncsSelection.imageMargin, model.basisFuncsSelection.noiseMargin);
              else
                dirs ~= "TXP%s_inr%s_%s_orth%s".format(model.pa.TX_POWER, model.INR, methodName, model.orthogonalizer.numOfTrainingSymbols);
            }

            foreach(i; 0 .. models.length)
                taskList.append((Model m, string dir){
                    JSONValue[] resList;

                    // 最初の一回は普通にやる
                    resList ~= mainImpl!methodName(m, dir);

                    static if(methodName.startsWith("SFHF"))
                    {
                        // writeln(mainImpl!methodName(m, dir)["training_symbols_per_second"]);
                        enum K = 10;    // 試行回数
                        uint sumOfSuccFreq;
                        JSONValue[] selectingRatioList;
                        foreach(j; 0 .. K){
                            m.rndSeed += 100;   // seed値を100ずつ足していく
                            auto res = mainImpl!methodName(m, null);
                            resList ~= res;
                            auto cnt = res["filterSpec"]["selectingIsSuccess"].array.map!(a => a.type == JSON_TYPE.TRUE ? 1 : 0).sum();
                            sumOfSuccFreq += cnt;
                            selectingRatioList ~= JSONValue(cnt / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling));
                        }

                        JSONValue jv = cast(JSONValue[string])null;

                        jv["cancellations"] = resList.map!(a => a["cancellation_dB"]).array();
                        
                        jv["selecting"] = (){
                            JSONValue[string] vv;
                            vv["selectingRatio"] = sumOfSuccFreq / cast(float)(m.ofdm.numOfFFT * m.ofdm.scaleOfUpSampling * K);
                            vv["selectingRatioList"] = selectingRatioList;

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
                            }
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

                            return vv;
                        }();
                        
                        auto file = File(buildPath(dir, "allResult.json"), "w");
                        file.write(jv.toPrettyString(JSONOptions.specialFloatLiterals));
                    }

                }, models[i], buildPath("results", dirs[i]));
            // foreach(i; 0 .. models.length)
            //     mainImpl!methodName(models[i], buildPath("results", dirs[i]));
        }

    //auto scheduler = new MPITaskScheduler();
    jobRun(taskList);
}


void main()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    // jobRun(1, 0, {
        mainJob();
    // });
}
