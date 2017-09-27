module simmain;

import std.algorithm;
import std.complex;
import std.datetime;
import std.json;
import std.math;
import std.mathspecial;
import std.numeric;
import std.path;
import std.range;
import std.random;
import std.stdio;

import carbon.math : nextPowOf2;

import dffdd.blockdiagram.amplifier;
import dffdd.blockdiagram.utils;
import dffdd.dsp.statistics;
import dffdd.utils.fft;
import dffdd.utils.unit;

import models;
import snippet;
import sigsim;


real qfunc(real x) { return 0.5 * erfc(x / SQRT2); }


auto psdSaveTo(R)(R r, string filename, string resultDir, size_t dropSize, Model model, bool* endFlag = null)
if(isForwardRange!R)
{
    alias E = ElementType!R;
    return r.loggingTo(makeSpectrumAnalyzer!E(filename, resultDir, dropSize, model, endFlag)).toWrappedRange;
}


auto makeSpectrumAnalyzer(C)(string filename, string resultDir, size_t dropSize, Model model, bool* endFlag = null)
{
    if(resultDir !is null) {
        return makeInstrument!C(delegate void(FiberRange!C r){
            r.drop(dropSize).writePSD(File(buildPath(resultDir, filename), "w"), model.samplingFreq, 1024);
            if(endFlag) *endFlag = true;
        });
    }else{
        if(endFlag) *endFlag = true;
        return makeInstrument!C(delegate void(FiberRange!C r){});
    }
}


auto makeCancellationProbe(C)(float* dst, string filename, string resultDir, size_t dropSize, Model model, bool* endFlag = null)
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

        if(endFlag) *endFlag = true;
    });
}


auto makeInBandCancellationProbe(C)(float* dst, string filename, string resultDir, size_t dropSize, Model model, bool* endFlag = null)
{
    alias Tup = Tuple!(C, C);

    return makeInstrument!Tup(delegate void(FiberRange!Tup r){
        auto rd = r.drop(dropSize);
        auto ratio = rd.calculateInBandSIC(model.samplingFreq, 1024, model.ofdm.numOfFFT, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling);
        if(dst !is null) *dst = 10*log10(ratio);

        if(resultDir !is null) {
            auto file = File(buildPath(resultDir, filename), "w");
            file.writefln("%s", 10*log10(ratio));
        }

        if(endFlag) *endFlag = true;
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

    with(model){
        useDesiredBaseband = false;
        useTxBaseband = true;           // 使う
        useTxBasebandSWP = true;        // 使う
        useDesiredPADirect = false;
        usePADirect = false;
        usePADirectSWP = false;
        useReceivedDesired = false;
        useReceivedSI = true;           // 使う
        useReceivedSISWP = true;        // 使う
        useReceived = true;             // 使う
        useNoise = true;                // 使う
    }

    auto signals = makeSimulatedSignals(model, resultDir);
    signals.trainAGC();

    // フィルタの学習
  enum string filterStructure = filterType.split("_")[0];
  enum string filterOptimizer = filterType.split("_")[1];

  enum bool isOrthogonalized = filterStructure[0] == 'O';

  static if(filterStructure.startsWith("PreIQI"))
  {
    import dffdd.filter.freqdomain;

    static if(filterStructure.startsWith("PreIQI-PH"))
        auto subfilter = makeParallelHammersteinFilter!(filterOptimizer, defaultDistortionOrder, false, false)(modOFDM(model), model);
    else static if(filterStructure.startsWith("PreIQI-FHF"))
        auto subfilter = makeFrequencyHammersteinFilter2!(filterOptimizer, 3, false)(model, false, false);
    else static assert(0);

    auto filter = new PreIQInvertionCanceller!(Complex!float, typeof(subfilter))(min(4, model.learningSymbols), model.ofdm.subCarrierMap, model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.scaleOfUpSampling, subfilter);
  }
  else static if(filterStructure.endsWith("PHDCM"))
    auto filter = makeParallelHammersteinWithDCMethodFilter!isOrthogonalized(modOFDM(model), model);
  else static if(filterStructure.endsWith("ARPH"))
    auto filter = makeAliasRemovableParallelHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("PH"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, defaultDistortionOrder, true, isOrthogonalized)(modOFDM(model), model);
  else static if(filterStructure.endsWith("CH"))
    auto filter = makeCascadeHammersteinFilter!(filterOptimizer)(modOFDM(model), model);
//   else static if(filterStructure.endsWith("CWLH"))
//     auto filter = makeCascadeWLHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
//   else static if(filterStructure.endsWith("CWL1H"))
//     auto filter = makeCascadeWL1HammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("IQISICFHF"))
  {
    import dffdd.filter.freqdomain;
    auto filter = new IQInversionSuccessiveInterferenceCanceller!(Complex!float, (defaultDistortionOrder+1)/2)(model.learningSymbols, 2, model.ofdm.subCarrierMap, model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.scaleOfUpSampling);
  }
  else static if(filterStructure.endsWith("WLFHF"))
  {
    static assert(!isOrthogonalized);
    auto filter = makeFrequencyHammersteinFilter2!(filterOptimizer, 1)(model, filterStructure.endsWith("SWLFHF"));
  }
  else static if(filterStructure.endsWith("DCMFHF"))
  {
    static assert(!isOrthogonalized);
    enum string filterOption = filterStructure[0 .. $-6];

    static assert(filterOption.canFind('1') || filterOption.canFind('2'));
    static assert(filterOption.canFind('P') || filterOption.canFind('C'));

    enum size_t type = filterOption.canFind('1') ? 1 : 2;
    enum Flag!"isParallel" isParallel = filterOption.canFind('P') ? Yes.isParallel : No.isParallel;

    static if(type == 2)
        auto filter = makeFrequencyDCMHammersteinFilter2!(type, isParallel, filterOptimizer)(model, filterStructure.endsWith("SDCMFHF"));
    else
        auto filter = makeFrequencyDCMHammersteinFilter!(type, isParallel, filterOptimizer)(model, filterStructure.endsWith("SDCMFHF"));
  }
  else static if(filterStructure.endsWith("CFHF"))
    static assert(0); // auto filter = makeFrequencyCascadeHammersteinFilter!(true, filterOptimizer)(model);
  else static if(filterStructure.endsWith("FHF")){
    static assert(!isOrthogonalized);
    auto filter = makeFrequencyHammersteinFilter2!(filterOptimizer)(model,
        filterStructure.endsWith("SFHF") || filterStructure.endsWith("S1FHF") || filterStructure.endsWith("S2FHF"),
        filterStructure.endsWith("S2FHF"));
  }else static if(filterStructure.endsWith("WL"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, 1)(modOFDM(model), model);
  else static if(filterStructure.endsWith("L"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, 1, false)(modOFDM(model), model);
  else static if(filterStructure.endsWith("TAYLOR"))
    auto filter = makeTaylorApproximationFilter!(1, false)(model);
  else
    static assert("Cannot identify filter model.");

    // フィルタの学習
    auto recvs = new Complex!float[model.blockSize],
         refrs = new Complex!float[model.blockSize],
         outps = new Complex!float[model.blockSize];

    {
        
        //assert(model.rndSeed != 893);

        filter.preLearning(model, (Model m){
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

        bool*[] endFlags = iota(5).map!"new bool"().array();
        auto psdBeforePSD = makeSpectrumAnalyzer!(Complex!float)("psd_beforeSIC.csv", resultDir, 0, model, endFlags[0]);
        auto psdAfterSIC = makeSpectrumAnalyzer!(Complex!float)("psd_afterSIC.csv", resultDir, 0, model, endFlags[1]);
        auto foutPSD = makeSpectrumAnalyzer!(Complex!float)("psd_filter_output.csv", resultDir, 0, model, endFlags[2]);
        auto sicValue = makeCancellationProbe!(Complex!float)(&sicv, "cancellation_value.csv", resultDir, 0, model, endFlags[3]);
        auto inbandSICValue = makeInBandCancellationProbe!(Complex!float)(null, "inband_cancellation_value.csv", resultDir, 0, model, endFlags[4]);

        StopWatch sw;
        size_t cancCNT;
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
                .put(inbandSICValue, recvs.zip(outps));
            }

            ++cancCNT;

            if(endFlags[].fold!"a && *b"(true))
                break;
        }

        infoResult["canceling_symbols_per_second"] = cancCNT * model.blockSize / model.ofdm.numOfSamplesOf1Symbol / ((cast(real)sw.peek.usecs) / 1_000_000);
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
        std.file.write(buildPath(resultDir, "info.json"), infoResult.toPrettyString(JSONOptions.specialFloatLiterals));

    return infoResult;
}
