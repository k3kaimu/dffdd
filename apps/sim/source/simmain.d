module simmain;

import std.algorithm;
import std.complex;
import core.time;
import std.datetime.stopwatch;
import std.json;
import std.math;
import std.mathspecial;
import std.numeric;
import std.path;
import std.range;
import std.random;
import std.stdio;
import std.typecons;
import std.conv;

import carbon.math : nextPowOf2;

import dffdd.blockdiagram.amplifier;
import dffdd.blockdiagram.utils;
import dffdd.dsp.statistics;
import dffdd.utils.fft;
import dffdd.utils.unit;
import dffdd.filter.primitives;
import dffdd.mod.primitives;

import models;
import snippet;
import sigsim;


double qfunc(double x) { return 0.5 * erfc(x / SQRT2); }


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


auto makeSyncSpectrumAnalyzer(C, size_t dim, alias names)(string filename, string resultDir, Model model, bool* endFlag = null)
{
    import std.meta : Repeat;
    alias Tup = Tuple!(Repeat!(dim, C));
    enum size_t N = 1024;

    if(resultDir !is null) {
        return makeInstrument!Tup(delegate void(FiberRange!Tup r){
            auto fftw = makeFFTWObject!C(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling);
            alias R = typeof(C.init.re);
            immutable size_t cp = model.ofdm.numOfCP * model.ofdm.scaleOfUpSampling;

            R[][] freq = new R[][](dim, fftw.inputs!R.length);
            foreach(buf; freq) foreach(ref e; buf) e = 0;

            foreach(ch; r.chunks(model.ofdm.numOfSamplesOf1Symbol).take(N)) {
                auto arr = ch.array();
                static foreach(i, name; names){{
                    arr[cp .. $].map!(a => a[i]).copy(fftw.inputs!R[]);
                    fftw.fft!R();
                    auto ps = fftw.outputs!R.map!(std.complex.sqAbs);
                    foreach(f, ref e; freq[i])
                        e += ps[f];
                }}
            }

            File file = File(buildPath(resultDir, filename), "w");
            file.writefln("Frequency,%(%s,%)", names);
            foreach(f; 0 .. fftw.inputs!R.length) {
                file.writef("%s,", f);
                foreach(i; 0 .. dim)
                    file.writef("%s,", freq[i][f] / N);
                file.writeln();
            }

            if(endFlag) *endFlag = true;
        });
    } else {
        if(endFlag) *endFlag = true;
        return makeInstrument!Tup(delegate void(FiberRange!Tup){});
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


auto makeFilter(string filterType)(Model model)
{
    enum string filterStructure = filterType.split("_")[0];
    enum string filterOptimizer = filterType.split("_")[1];

    enum bool isOrthogonalized = filterStructure[0] == 'O';

  static if(filterStructure.startsWith("PreIQI"))
  {
    import dffdd.filter.freqdomain;

    static if(filterStructure.startsWith("PreIQI-PH"))
        auto subfilter = makeParallelHammersteinFilter!(filterOptimizer, defaultDistortionOrder, false, false)(modOFDM(model), model);
    else static if(filterStructure.startsWith("PreIQI-FHF"))
        auto subfilter = makeFrequencyHammersteinFilter2!(filterOptimizer, 3, false)(model);
    else static if(filterStructure.startsWith("PreIQI-L"))
        auto subfilter = makeParallelHammersteinFilter!(filterOptimizer, 1, false, false)(modOFDM(model), model);
    else static assert(0);

    auto filter = new PreIQInvertionCanceller!(Complex!float, typeof(subfilter))(min(4, model.learningSymbols), model.ofdm.subCarrierMap, model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.scaleOfUpSampling, subfilter);
  }
  else static if(filterStructure.endsWith("PHDCM"))
    auto filter = makeParallelHammersteinWithDCMethodFilter!isOrthogonalized(modOFDM(model), model);
  else static if(filterStructure.endsWith("ARPH"))
    auto filter = makeAliasRemovableParallelHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure[0 .. $-1].endsWith("SubPH"))
  {
    enum size_t POrder = filterStructure[$-1 .. $].to!int;
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, SubSetOfPADistorter!(Complex!float, POrder), isOrthogonalized)(modOFDM(model), model);
  }
  else static if(filterStructure[0 .. $-1].endsWith("PHPAOnly"))
  {
    enum size_t POrder = filterStructure[$-1 .. $].to!int;
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, OnlyPADistorter!(Complex!float, POrder), isOrthogonalized)(modOFDM(model), model);
  }
  else static if(filterStructure[0 .. $-1].endsWith("PH"))
  {
    enum size_t POrder = filterStructure[$-1 .. $].to!int;
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, PADistorter!(Complex!float, POrder), isOrthogonalized)(modOFDM(model), model);
  }
  else static if(filterStructure.endsWith("CH"))
    auto filter = makeCascadeHammersteinFilter!(filterOptimizer)(modOFDM(model), model);
  // else static if(filterStructure.endsWith("CWLH"))
  //   auto filter = makeCascadeWLHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  // else static if(filterStructure.endsWith("CWL1H"))
  //   auto filter = makeCascadeWL1HammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("IterativeFreqSIC"))
  {
    import dffdd.filter.freqdomain;
    auto filter = new IQInversionSuccessiveInterferenceCanceller!(Complex!float, (defaultDistortionOrder+1)/2)(model.learningSymbols, model.iterativeFreqSIC.iterations, model.ofdm.subCarrierMap, model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.scaleOfUpSampling);
  }
  else static if(filterStructure[0 .. $-1].endsWith("Sidelobe"))
  {
    enum size_t POrder = filterStructure[$-1 .. $].to!int / 2 + 1;

    import dffdd.filter.sidelobe;
    auto filter = new SidelobeIterativeWLNL!(Complex!float, POrder)(
        model.learningSymbols, model.iterativeFreqSIC.iterations,
        model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling,
        model.channel.taps, No.isChFreqEst, Yes.isInvertRX,
        Yes.useNewton, model.iterativeFreqSIC.newtonIterations, model.iterativeFreqSIC.use3rdSidelobe,
        model.iterativeFreqSIC.numOfSCForEstNL, model.iterativeFreqSIC.estimationOrder);
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
  else static if(filterStructure[0 .. $-1].endsWith("FHF"))
  {
    // static assert(!isOrthogonalized);

    enum size_t POrder = filterStructure[$-1 .. $].to!int;

    static if(filterStructure.canFind("Sub"))
        alias Dist = SubSetOfPADistorter!(Complex!float, POrder);
    else static if(filterStructure.canFind("PHPAOnly"))
        alias Dist = OnlyPADistorter!(Complex!float, POrder);
    else
        alias Dist = PADistorter!(Complex!float, POrder);

    static if(isOrthogonalized)
    {
        import dffdd.filter.orthogonalize;

        alias GS = GramSchmidtOBFFactory!(Complex!float);
        auto dist = new OrthogonalizedVectorDistorter!(Complex!float, Dist, GS)(new Dist(), new GS(Dist.outputDim));
    }
    else
    {
        auto dist = new Dist();
    }

    auto freqFilter = makeFrequencyHammersteinFilter2!(filterOptimizer, typeof(dist))(dist, model);

    static if(filterStructure.canFind("S2"))
    {
        auto filter = makeFrequencyDomainBasisFunctionSelector(model, freqFilter);
    }
    else
    {
        alias filter = freqFilter;
    }
  }else static if(filterStructure.endsWith("WL"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, PADistorter!(Complex!float, 1), isOrthogonalized)(modOFDM(model), model);
  else static if(filterStructure.endsWith("L"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, OnlyPADistorter!(Complex!float, 1), false)(modOFDM(model), model);
  else static if(filterStructure.endsWith("TAYLOR"))
    auto filter = makeTaylorApproximationFilter!(1, false)(model);
  else static if(filterStructure.endsWith("Nop"))
    auto filter = new NopCanceller!(Complex!float)();
  else static if(filterStructure.endsWith("Ahmed2013"))
  {
    import dffdd.filter.iterative_ahmed_2013;
    auto filter = new IterativeAhmed2013!(Complex!float)(
        model.learningSymbols, model.iterativeFreqSIC.iterations, model.iterativeFreqSIC.newtonIterations,
        model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling,
        model.channel.taps
    );
  }
  else
    static assert("Cannot identify filter model.");

    pragma(msg, typeof(filter).stringof ~ " is instantiated.");

    return filter;
}


JSONValue mainImpl(string filterType)(Model model, string resultDir = null)
{
    enum string filterStructure = filterType.split("_")[0];
    enum string filterOptimizer = filterType.split("_")[1];


    JSONValue infoResult = ["type": filterType];

    immutable bool* alwaysFalsePointer = new bool(false);
    immutable bool* alwaysTruePointer = new bool(true);

    auto fftObject = makeFFTWObject!Complex(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling);

    immutable ofdmModSignalPower = (){
        auto _modOFDMTest = modOFDM(model);
        return randomBits(1, model).connectToModulator(_modOFDMTest, alwaysFalsePointer, model).measurePower(1024*1024);
    }();

    if(resultDir !is null){
        import std.file : mkdirRecurse;
        mkdirRecurse(resultDir);
    }

    auto signals = makeSimulatedSignals(model, resultDir);
    signals.trainAGC();

    // フィルタを作成．信号のLog出力をするのであれば，ラッパをかませる．
    static if(filterType.startsWith("Log_"))
    {
        import std.file : mkdirRecurse;
        mkdirRecurse("signallog");
        auto filter = makeSignalLoggingCanceller!(Complex!float)(makeFilter!(filterType[4 .. $])(model), "signallog", model.uniqueId);
    }
    else
        auto filter = makeFilter!filterType(model);

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

    signals.ignoreDesired = true;
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
            static if(filterStructure[0 .. $-1].endsWith("FHF") || filterStructure.endsWith("IterativeFreqSIC"))
            {
                if(blockIdx >= model.swappedSymbols * model.ofdm.numOfSamplesOf1Symbol / model.blockSize)
                    signals.useSWPOFDM = false;
                else
                    signals.useSWPOFDM = true;
            }

            signals.fillBuffer!(["txBaseband", "received"])(refrs, recvs);

            sw.start();
            filter.apply!(Yes.learning)(refrs, recvs, outps);
            sw.stop();

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


        infoResult["training_symbols_per_second"] = model.numOfFilterTrainingSymbols / ((cast(real)sw.peek.total!"usecs") / 1_000_000);
    }

    {
        // float sicv;

        bool*[] endFlags = iota(5).map!"new bool"().array();
        auto psdBeforePSD = makeSpectrumAnalyzer!(Complex!float)("psd_beforeSIC.csv", resultDir, 0, model, endFlags[0]);
        auto psdAfterSIC = makeSpectrumAnalyzer!(Complex!float)("psd_afterSIC.csv", resultDir, 0, model, endFlags[1]);
        auto foutPSD = makeSpectrumAnalyzer!(Complex!float)("psd_filter_output.csv", resultDir, 0, model, endFlags[2]);
        // auto sicValue = makeCancellationProbe!(Complex!float)(&sicv, "cancellation_value.csv", resultDir, 0, model, endFlags[3]);
        auto inbandSICValue = makeInBandCancellationProbe!(Complex!float)(null, "inband_cancellation_value.csv", resultDir, 0, model, endFlags[3]);
        auto syncSpecAnalyzer = makeSyncSpectrumAnalyzer!(Complex!float, 2, ["Y", "D"])("psd_sync.csv", resultDir, model, endFlags[4]);

        double sumRemainPower = 0;
        double sumSI = 0;
        StopWatch sw;
        size_t cancCNT;
        foreach(blockIdxo; 0 .. 1024)
        {
            signals.fillBuffer!(["txBaseband", "received"])(refrs, recvs);

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
                // .put(sicValue, recvs.zip(outps));
                .put(inbandSICValue, recvs.zip(outps));
                .put(syncSpecAnalyzer, zip(recvs, outps));
            }

            sumRemainPower += outps.map!(a => a.sqAbs).sum();
            sumSI += recvs.map!(a => a.sqAbs).sum();
            ++cancCNT;

            if(endFlags[].fold!"a && *b"(true))
                break;
        }

        sumRemainPower /= cancCNT * model.blockSize;
        sumSI /= cancCNT * model.blockSize;
        infoResult["canceling_symbols_per_second"] = cancCNT * model.blockSize / model.ofdm.numOfSamplesOf1Symbol / ((cast(double)sw.peek.total!"usecs") / 1_000_000);
        // infoResult["cancellation_dB"] = sicv;

        // 雑音電力を測る
        signals.ignoreSI = true;
        scope(exit) signals.ignoreSI = false;
        double sumNoisePower = 0;
        foreach(i; 0 .. cancCNT) {
            signals.fillBuffer!(["received"])(refrs);
            sumNoisePower += refrs.map!(a => a.sqAbs).sum();
        }
        sumNoisePower /= cancCNT * model.blockSize;
        double sumRemainSI = sumRemainPower - sumNoisePower;
        if(sumRemainSI < 0) sumRemainSI = 0;
        Gain RINR = Gain.fromPowerGain(sumRemainSI / sumNoisePower),
             INR = Gain.fromPowerGain(sumSI / sumNoisePower),
             canc = Gain.fromPowerGain((INR.asP + 1) / (RINR.asP + 1));

        infoResult["RemainPower"] = sumRemainPower;
        infoResult["SIPower"] = sumSI;
        infoResult["cancellation_dB"] = canc.asdB;
        infoResult["RINR_dB"] = RINR.asdB;
        infoResult["INR_dB"] = INR.asdB;
        if(resultDir !is null){
            File(buildPath(resultDir ,"INR_value.csv"), "w").writeln(INR.asdB);
            File(buildPath(resultDir ,"cancellation_value.csv"), "w").writeln(canc.asdB);
            File(buildPath(resultDir ,"RINR_value.csv"), "w").writeln(RINR.asdB);
        }

    }

    if(model.outputBER)
        simulateMeasureBEREVMImpl(filter, signals, model, infoResult, resultDir, fftObject);

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

    infoResult["modelSpec"] = signals.info();

    if(resultDir !is null){
        import std.file : filewrite = write;
        filewrite(buildPath(resultDir, "info.json"), infoResult.toPrettyString(JSONOptions.specialFloatLiterals));
    }

    return infoResult;
}



void simulateMeasureBEREVMImpl(Filter, Signals, FFTObject)(ref Filter filter, ref Signals signals, ref Model model, ref JSONValue infoResult, string resultDir, ref FFTObject fftObject)
{
    import dffdd.mod.qam;
    import dffdd.mod.ofdm;
    import dffdd.mod.primitives : Bit;

    alias F = float;
    alias C = Complex!float;

    immutable nCP = model.ofdm.numOfCP * model.ofdm.scaleOfUpSampling;
    immutable nFFT = model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling;
    immutable nSC = model.ofdm.numOfSubcarrier;
    immutable nSYM = model.ofdm.numOfSamplesOf1Symbol;

    auto ofdmMod = new dffdd.mod.ofdm.OFDM!(Complex!float)(model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling);
    auto qamMod = QAM!(Complex!float)(16);

    if(!model.withSI) signals.ignoreSI = true;
    signals.ignoreDesired = false;
    scope(exit){
        signals.ignoreSI = false;
        signals.ignoreDesired = true;
    }

    // // model.channelDesiredからチャネルの周波数応答の真値を得る．
    // auto channelFreqResp = (){
    //     C[] buf = new C[nSYM];
    //     C[] dst = new C[nSC];
    //     buf[] = C(0);
    //     buf[nCP .. nCP + model.channelDesired.taps] = model.channelDesired.impulseResponse.map!(a => C(a)).array()[];
    //     ofdmMod.demodulate(buf, dst);

    //     auto g = signals.desiredSignalGain.asV;
    //     real scale = sqrt((nFFT * 1.0)^^2 / nSC);   // ofdmMod.demodulateの中でかかる係数
    //     foreach(ref e; dst)
    //         e *= g * scale;

    //     return dst;
    // }();
    C[] channelFreqResp;
    C[] channelFreqCorr, txFreqSqAbs;
    channelFreqCorr = new C[](nSC);
    txFreqSqAbs = new C[](nSC);
    channelFreqCorr[] = C(0);
    txFreqSqAbs[] = C(0);
    size_t numDesiredChannelTraining;
    immutable maxNumDesiredChannelTraining = 10;    // 何シンボルでチャネルを学習するか

    C[] receivedSCs = new C[nSC],
        referenceSCs = new C[nSC];

    void trainingDesiredChannelFreqResp(in C[] receivedSymbol, in C[] referenceSymbol)
    {
        if(channelFreqResp !is null) return;

        ofdmMod.demodulate(receivedSymbol, receivedSCs);
        ofdmMod.demodulate(referenceSymbol, referenceSCs);

        foreach(i; 0 .. receivedSCs.length) {
            channelFreqCorr[i] += receivedSCs[i] * referenceSCs[i].conj;
            txFreqSqAbs[i] += referenceSCs[i].sqAbs;
        }

        ++numDesiredChannelTraining;

        if(numDesiredChannelTraining >= maxNumDesiredChannelTraining) {
            channelFreqResp = new C[](nSC);
            foreach(i, ref e; channelFreqResp)
                e = channelFreqCorr[i] / txFreqSqAbs[i];
        }
    }

    size_t totalOFDMSymbols;

    ushort[] receivedQAMSyms, referenceQAMSyms;
    BERCounter counter = BERCounter(qamMod.symInputLength);

    auto sumPs = new double[nSC],
         diffPs = new double[nSC];
    sumPs[] = 0;
    diffPs[] = 0;

    File receivedSCOutput, referenceSCOutput;
    if(resultDir !is null) {
        receivedSCOutput = File(buildPath(resultDir, "receivedSymbols.csv"), "w");
        referenceSCOutput = File(buildPath(resultDir, "referenceSymbols.csv"), "w");
    }


    void checkBERAndEVMForEachSymbol(in C[] receivedSymbol, in C[] referenceSymbol)
    {
        ofdmMod.demodulate(receivedSymbol, receivedSCs);
        ofdmMod.demodulate(referenceSymbol, referenceSCs);

        // 通ってきたチャネルを等化する
        foreach(i; 0 .. nSC)
            receivedSCs[i] /= channelFreqResp[i];

        if(totalOFDMSymbols < 300 && resultDir !is null) {
            receivedSCOutput.writefln("%(%s,%)", receivedSCs.map!"[a.re,a.im]".joiner);
            referenceSCOutput.writefln("%(%s,%)", referenceSCs.map!"[a.re,a.im]".joiner);
        }

        // QAMの復調をする
        qamMod.demodulate_symbol(receivedSCs, receivedQAMSyms);
        qamMod.demodulate_symbol(referenceSCs, referenceQAMSyms);

        counter.count(receivedQAMSyms, referenceQAMSyms);

        foreach(i; 0 .. nSC){
            diffPs[i] += (receivedSCs[i] - referenceSCs[i]).sqAbs;
            sumPs[i] += referenceSCs[i].sqAbs;
        }

        totalOFDMSymbols += 1;
    }

    auto rcvAfterSICPSD = makeSpectrumAnalyzer!(Complex!float)("psd_rcv_afterSIC.csv", resultDir, 0, model);
    auto rcvBeforeSICPSD = makeSpectrumAnalyzer!(Complex!float)("psd_rcv_beforeIC.csv", resultDir, 0, model);
    
    auto receivedSymbol = new C[nSYM],
         txSymbol = new C[nSYM],
         canceledSymbol = new C[nSYM],
         desiredReferenceSymbol = new C[nSYM];

    size_t loopcount = 0;
    while(1){
        signals.fillBuffer!(["txBaseband", "received", "desiredBaseband"])(txSymbol, receivedSymbol, desiredReferenceSymbol);

        if(model.withSIC && model.withSI)
            filter.apply!(No.learning)(txSymbol, receivedSymbol, canceledSymbol);
        else
            canceledSymbol[] = receivedSymbol[];

        // 最初の10シンボルは飛ばす
        ++loopcount;
        if(loopcount <= 10)
            continue;

        if(channelFreqResp is null)
            trainingDesiredChannelFreqResp(canceledSymbol, desiredReferenceSymbol);
        else
            checkBERAndEVMForEachSymbol(canceledSymbol, desiredReferenceSymbol);

        if(counter.totalBits >= model.berCounter.totalBits
        && totalOFDMSymbols >= model.berCounter.evmSymbols)
            break;
    }

    with(counter.result) {
        infoResult["ber"] = ber;
        infoResult["ser"] = ser;
    }
    infoResult["evm"] = diffPs.zip(sumPs).map!"a[0] / a[1]".array();
}
