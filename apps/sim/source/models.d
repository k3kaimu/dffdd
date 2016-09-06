module models;

import core.thread;

import std.algorithm;
import std.complex;
import std.conv;
import std.datetime;
import std.file : mkdirRecurse;
import std.format;
import std.functional;
import std.math;
import std.mathspecial;
import std.meta;
import std.parallelism;
import std.path;
import std.random;
import std.range;
import std.stdio;
import std.typecons;

import dffdd.blockdiagram.adder;
import dffdd.blockdiagram.amplifier;
import dffdd.blockdiagram.decimator;
import dffdd.blockdiagram.filter;
import dffdd.blockdiagram.iqmixer;
import dffdd.blockdiagram.mod.ofdm;
import dffdd.blockdiagram.noise;
import dffdd.blockdiagram.quantizer;
import dffdd.blockdiagram.txchain;
import dffdd.blockdiagram.utils;
import dffdd.filter.diagonal;
import dffdd.filter.lms;
import dffdd.filter.ls;
import dffdd.filter.mempoly;
import dffdd.filter.orthfreqpolyfil;
import dffdd.filter.orthogonalize;
import dffdd.filter.polynomial;
import dffdd.filter.polynomial;
import dffdd.filter.primitives;
import dffdd.filter.rls;
import dffdd.filter.state;
import dffdd.utils.fft;
import dffdd.utils.unit;
//import dffdd.utils.msgpackrpc;




//import carbon.range;
import dranges.range;
import dranges.algorithm;

import dffdd.itpp;
import dffdd.mod.primitives;
import dffdd.mod.bpsk;
import dffdd.mod.qpsk;
import dffdd.mod.qam;
import dffdd.mod.ofdm;

import dffdd.dsp.statistics;


import constant;
import snippet;

alias BasisFunctions = AliasSeq!(x => x,
                              x => x.conj,
                              x => x * (x.re^^2 + x.im^^2),
                              x => x.conj * (x.re^^2 + x.im^^2),
                              x => x * (x.re^^2 + x.im^^2)^^2,
                              x => x.conj * (x.re^^2 + x.im^^2)^^2,
                              //x => x * (x.re^^2 + x.im^^2)^^3,
                              //x => x.conj * (x.re^^2 + x.im^^2)^^3,
                              );

struct Model
{
    size_t numOfModelTrainingSymbols = 1024;
    size_t numOfFilterTrainingSymbols = 1024*10;
    //size_t blockSize = 1024;
    size_t blockSize() const @property { return ofdm.numOfSamplesOf1Symbol*4; }
    real samplingFreq = 20e6 * 4;
    real SNR = 20;
    real INR = 60;
    bool withSIC = true;
    bool withSI = true;
    size_t numOfPopFront = 1;
    size_t learningSymbols = 10;
    size_t learningCount = 10;


    bool useDTXIQ = true;
    bool useDTXPA = true;
    bool useDTXIQ2 = false;
    bool useSTXIQ = true;
    bool useSTXPA = true;
    bool useSTXIQ2 = false;
    bool useSRXLN = true;
    bool useSRXIQ = true;
    bool useSRXQZ = true;

/*
    struct QAM
    {
        uint arity = 16;
    }
    QAM qam;
*/

    struct OFDM
    {
        uint numOfFFT = 64;
        uint numOfCP = 16;
        uint numOfSubcarrier = 52;
        uint scaleOfUpSampling = 4;
        real PAPR = 10;                 // 10dB
        //uint model.numOfSamplesOf1Symbol = 
        uint numOfSamplesOf1Symbol() const @property { return scaleOfUpSampling * (numOfFFT + numOfCP); }
        bool[] subCarrierMap() const @property
        {
            bool[] ret = new bool[scaleOfUpSampling * numOfFFT];

            assert(numOfSubcarrier % 2 == 0);

            ret[1 .. numOfSubcarrier/2 + 1] = true;
            ret[$ - numOfSubcarrier/2 .. $] = true;

            return ret;
        }
    }
    OFDM ofdm;


    struct ThermalNoise
    {
        real temperature = 300;
        real power(Model model)
        {
            return noisePower(model.samplingFreq / (model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling) * model.ofdm.numOfSubcarrier,
                              model.thermalNoise.temperature);
        }
    }
    ThermalNoise thermalNoise;


    //struct UpDownSampler
    //{
    //    uint scaleOfUpSampling = 8;
    //    immutable(real)[] decimationFIRFilterTaps = [
    //        0.0064217,0.0030971,0.003501,0.0036559,0.0034831,0.0029178,0.0019183,0.00047411,-0.0013877,-0.003598,-0.0060443,-0.0085724,
    //        -0.010993,-0.013088,-0.014628,-0.015381,-0.015132,-0.0137,-0.010952,-0.0068179,-0.0013004,0.0055175,0.013472,0.022324,
    //        0.031764,0.041433,0.050934,0.059859,0.06781,0.074421,0.079383,0.082458,0.0835,0.082458,0.079383,0.074421,0.06781,
    //        0.059859,0.050934,0.041433,0.031764,0.022324,0.013472,0.0055175,-0.0013004,-0.0068179,-0.010952,-0.0137,-0.015132,
    //        -0.015381,-0.014628,-0.013088,-0.010993,-0.0085724,-0.0060443,-0.003598,-0.0013877,0.00047411,0.0019183,0.0029178,
    //        0.0034831,0.0036559,0.003501,0.0030971,0.0064217];

    //    size_t numOfShift = 67;
    //    //size_t numOfShift = 0;
    //}
    //UpDownSampler upDownSampler;


    struct PowerSpectralDensityAnalyzer
    {
        uint resolution = 1024;
    }
    PowerSpectralDensityAnalyzer powerSpectralDensityAnalyzer;


    struct BERCounter
    {
        ulong totalBits = 10_000_000;
    }
    BERCounter berCounter;


    struct TXIQMixer
    {
        real IIR = 25;  // 25dB
    }
    TXIQMixer txIQMixer;


    struct RXIQMixer
    {
        real IIR = 25;  // 25dB
    }
    RXIQMixer rxIQMixer;


    struct PA 
    {
        real GAIN = 27;     // 30dB
        real IIP3 = 17;     // 20dBm
        real TX_POWER = 27; // 20dBm
    }
    PA pa;

    struct LNA
    {
        real GAIN = 20;     // 20dB
        real NF = 4;        // 4dB
    }
    LNA lna;


    struct Quantizer
    {
        uint numOfBits = 24;
    }
    Quantizer quantizer;


    struct AWGN
    {
        real SNR = 24;      // 24dB
    }
    AWGN awgn;
}


auto randomBits(uint prn, Model)
{
    import dffdd.gps.code;
    return L2CLCode(prn).map!"cast(ubyte)(a < 0 ? 1 : 0)";
}


auto siRandomBits(Model model) { return randomBits(1, model); }
auto desiredRandomBits(Model model) { return randomBits(193, model); }



auto modOFDM(Model model)
{
    return chainedMod(
        //dffdd.mod.bpsk.BPSK.init,
        dffdd.mod.qam.QAM(16),
        new dffdd.mod.ofdm.OFDM!(Complex!float)(model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling),
    );
}


auto connectToModulator(R, Mod)(R r, Mod modObj, const(bool)* swExchange, Model model)
{
    auto swModObj = modObj.makeOFDMSymbolExchanger(swExchange, model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling);

    return r
    .splitN(swModObj.symInputLength)
    .tmap!(reverseArgs!mod, [0])(swModObj)
    .joiner()
    .toWrappedRange
    ;
}


auto connectToDemodulator(R, Mod)(R r, Mod modObj/*, Bits bits*/, Model)
{
    static
    Complex!float makeCpx(F)(Complex!F r) { return Complex!float(r.re, r.im); } 


    return r
    .map!makeCpx
    .splitN(modObj.symOutputLength)
    //.connectToOFDMEqualizer(modObj, bits)
    .tmap!(reverseArgs!demod, [0])(modObj)
    .joiner()
    .toWrappedRange
    ;
}


//auto connectToUpSampler(R)(R r, Model model)
//{
//    return r
//    .tmap!"a.repeat(b).array()"(constant.upDownSampler.scaleOfUpSampling).joiner()
//    .connectTo!FIRFilter(constant.upDownSampler.decimationFIRFilterTaps)
//    .toWrappedRange
//    ;
//}



//auto connectToDownSampler(R)(R r, Model model)
//{
//    return r
//    .drop(constant.upDownSampler.numOfShift)
//    .connectTo!FIRFilter(constant.upDownSampler.decimationFIRFilterTaps)
//    .connectTo!(SimpleDecimator!(constant.upDownSampler.scaleOfUpSampling))
//    .toWrappedRange
//    ;
//}


auto thermalNoise(Model model)
{
    return ThermalNoise(model.samplingFreq, model.thermalNoise.temperature);
}


auto connectToAWGN(R)(R r, Model model)
{
    return r.add(thermalNoise(model)).toWrappedRange;
}


//auto connectToFlatRayleighFadingChannel(R)(R r)
//{
//    BoxMuller!Random bm;
//    bm.seed(unpredictableSeed());
//    bm.popFront();

//    return r
//    .zip(
//        bm
//        .tmap!"(a/SQRT2).repeat(b)"((Constant.OFDM.numOfFFT + Constant.OFDM.numOfCP)*Constant.OFDM.scaleOfUpSampling*Constant.UpDownSampler.scaleOfUpSampling*100)
//        .joiner
//    )
//    .map!"a[0]*a[1]";
//}


//auto connectToOFDMEqualizer(R, Mod, Bits)(R r, Mod modObj, Bits bits)
//{
//    return OFDMEqualizer.makeBlock(
//        r,
//        modObj,
//        Constant.OFDM.numOfFFT,
//        Constant.OFDM.numOfCP,
//        Constant.OFDM.numOfSubcarrier,
//        Constant.OFDM.scaleOfUpSampling,
//        bits).toWrappedRange;
//}


auto connectToTXIQMixer(R)(R r, Model model)
{
    return r.connectTo!IQImbalance(0.dB, model.txIQMixer.IIR.dB).toWrappedRange;
}


auto connectToRXIQMixer(R)(R r, Model model)
{
    return r.connectTo!IQImbalance(0.dB, model.rxIQMixer.IIR.dB).toWrappedRange;
}


auto connectToPowerAmplifier(R)(R r, Model model)
{
    auto v = (model.pa.TX_POWER - model.pa.GAIN).dBm;

    return r
    .connectTo!PowerControlAmplifier(v)
    .connectTo!PowerAmplifier(model.pa.GAIN.dB, model.pa.IIP3.dBm)
    .toWrappedRange;
}


shared immutable(Complex!float)[] channelCoefsOfSI;

shared static this()
{
    Random rGen;
    rGen.seed(114514);

    Complex!float[] coefs;

    coefs ~= Complex!float(1, 1);

    foreach(i; 1 .. 8)
        coefs ~= Complex!float(uniform01(rGen), uniform01(rGen));

    foreach(i; 8 .. 16)
        coefs ~= Complex!float(uniform01(rGen), uniform01(rGen)) / sqrt(10.0f);

    foreach(i; 16 .. 32)
        coefs ~= Complex!float(uniform01(rGen), uniform01(rGen)) / sqrt(20.0f);

    foreach(i; 32 .. 64)
        coefs ~= Complex!float(uniform01(rGen), uniform01(rGen)) / sqrt(30.0f);

    channelCoefsOfSI = coefs.dup;
}


auto connectToMultiPathChannel(R)(R r)
{
    return r.connectTo!FIRFilter(channelCoefsOfSI).toWrappedRange;
}


auto connectToLNA(R)(R r, Model model)
{
    return r
    .add(thermalNoise(model).tmap!"a*b"(model.lna.NF.dB.gain - 1))
    .tmap!"a*b"(model.lna.GAIN.dB.gain)
    .toWrappedRange;
}


auto connectToQuantizer(R)(R r, Model model)
{
    return r
    .connectTo!PowerControlAmplifier((30 - model.ofdm.PAPR).dBm)
    .connectTo!SimpleQuantizer(model.quantizer.numOfBits)
    .toWrappedRange;
}


auto connectToTxChain(R)(R r, Model model)
{
    //return r.connectToTXIQMixer().connectToPowerAmplifier();
    return TXChain.makeBlock(r, 1/model.txIQMixer.IIR.dB.gain, 0, model.pa.GAIN.dB, model.pa.IIP3.dBm).toWrappedRange();
}


auto connectToRxChain(R)(R r)
{
    return r.connectToLNA().connectToRXIQMixer().connectToQuantizer();
}


auto makeParallelHammersteinFilter(bool isOrthogonalized, size_t numOfBasisFuncs = BasisFunctions.length, Mod)(Mod mod, Model model)
{
    alias BFs = BasisFunctions[0 .. numOfBasisFuncs];

    import dffdd.filter.diagonal;
    import dffdd.filter.lms;
    import dffdd.filter.ls;
    import dffdd.filter.rls;
    import dffdd.filter.mempoly;
    import dffdd.filter.polynomial;
    import dffdd.filter.orthogonalize;
    import dffdd.filter.state;

    //auto state = new MemoryPolynomialState!(Complex!float, 8, 4, 0, 0, false, true)(1);

    //writeln("Use");
    //auto adapter = new LMSAdapter!(typeof(state))(state, 0.001, 1024, 0.5);
    //auto adapter = makeRLSAdapter(state, 1 - 1E-4, 1E-7);
    //auto adapter = lsAdapter(state, 10000);

    //return new PolynomialFilter!(typeof(state), typeof(adapter))(state, adapter);

    auto orthogonalizer = new GramSchmidtOBFFactory!(Complex!float, BFs)();
    //auto orthogonalizer = new DiagonalizationOBFFactory!(Complex!float, BFs)();

    {
        orthogonalizer.start();
        scope(exit)
            orthogonalizer.finish();


        auto signal = randomBits(1, model).connectToModulator(mod, new bool(false), model);

        //.put(orthogonalizer, signal.take(1024*400));
        Complex!float[] buf = new Complex!float[1024];
        foreach(i; 0 .. 1024){
            foreach(j; 0 .. 1024){
                auto f = signal.front;
                buf[j] = complex(f.re, f.im);
                signal.popFront();
            }

            .put(orthogonalizer, buf);
        }
    }

    Complex!float[][BFs.length] coefs;
    foreach(i, ref e; coefs){
        e = new Complex!float[BFs.length];
        orthogonalizer.getCoefs(i, e);
    }


    Complex!float delegate(Complex!float x)[BFs.length] bflist;
    foreach(i, BF; BFs){
      static if(isOrthogonalized)
        bflist[i] = delegate Complex!float (Complex!float x) { return OBFEval!BFs(x, coefs[i]); };
      else
        bflist[i] = delegate Complex!float (Complex!float x) { return BF(x); };
    }
        //bflist[i] = delegate Complex!float (Complex!float x) { return BF(x); };

    auto state = new ParallelHammersteinState!(Complex!float, BFs.length, true)(64, bflist);
    //auto adapter = new LMSAdapter!(typeof(state))(state, 0.001, 1024, 0.5);
    //auto adapter = makeRLSAdapter(state, 1 - 1E-4, 1E-7);
    immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;
    auto adapter = lsAdapter(state, 80 * 4 * model.learningSymbols).trainingLimit(samplesOfOnePeriod * model.learningCount);

    return polynomialFilter(state, adapter);
}


/**
Fast Training and Cancelling method for Hammerstein Self-Interference Canceller on Full-Duplex OFDM Communication
*/
auto makeCascadeHammersteinFilter(bool isOrthogonalized, alias filterBuilder = serialFilter, Mod)(Mod mod, Model model)
{
    import dffdd.filter.diagonal;
    import dffdd.filter.lms;
    import dffdd.filter.ls;
    import dffdd.filter.rls;
    import dffdd.filter.mempoly;
    import dffdd.filter.polynomial;
    import dffdd.filter.orthogonalize;
    import dffdd.filter.state;
    import std.meta;
    import std.stdio;

    //writeln("orthogonalizer start");

    auto orthogonalizer = new GramSchmidtOBFFactory!(Complex!float, BasisFunctions)();
    //auto orthogonalizer = new DiagonalizationOBFFactory!(Complex!float, BasisFunctions)();

    {
        orthogonalizer.start();
        scope(exit)
            orthogonalizer.finish();

        auto signal = randomBits(1, model).connectToModulator(mod, new bool(false), model);

        //.put(orthogonalizer, signal.take(1024*400));
        Complex!float[] buf = new Complex!float[1024];
        foreach(i; 0 .. 1024){
            foreach(j; 0 .. 1024){
                auto f = signal.front;
                buf[j] = complex(f.re, f.im);
                signal.popFront();
            }
            .put(orthogonalizer, buf);
        }
    }


  static if(isOrthogonalized)
  {
    pragma(msg, "Orthogonalized");
    Complex!float[][BasisFunctions.length] coefs;
    foreach(i, ref e; coefs){
        e = new Complex!float[BasisFunctions.length];
        orthogonalizer.getCoefs(i, e);
    }


    auto adaptInputTransformer(size_t i, S)(S state)
    {
        return state.inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[i].dup);
    }
  }
  else
  {
    pragma(msg, "Un-Orthogonalized");
    auto adaptInputTransformer(size_t i, S)(S state)
    {
        return state.inputTransformer!(x => BasisFunctions[i](x))();
    }
  }

    auto st1 = adaptInputTransformer!0(new FIRState!(Complex!float, true)(64));
    auto st12 = adaptInputTransformer!0(new FIRState!(Complex!float, true)(2));
    auto st1c = adaptInputTransformer!1(new FIRState!(Complex!float, true)(64));
    auto st12c = adaptInputTransformer!1(new FIRState!(Complex!float, true)(2));
    //auto st2 = adaptInputTransformer!2(new FIRState!(Complex!float, true)(64));
    auto st3 = adaptInputTransformer!2(new FIRState!(Complex!float, true)(64));
    auto st32 = adaptInputTransformer!2(new FIRState!(Complex!float, true)(2));
    //auto st32 = adaptInputTransformer!2(new FIRState!(Complex!float, true)(2));
    auto st3c = adaptInputTransformer!3(new FIRState!(Complex!float, true)(64));
    auto st32c = adaptInputTransformer!3(new FIRState!(Complex!float, true)(2));
    auto st5 = adaptInputTransformer!4(new FIRState!(Complex!float, true)(64));
    auto st52 = adaptInputTransformer!4(new FIRState!(Complex!float, true)(2));
    auto st5c = adaptInputTransformer!5(new FIRState!(Complex!float, true)(64));
    auto st52c = adaptInputTransformer!5(new FIRState!(Complex!float, true)(2));
    //auto st7 = new FIRState!(Complex!float, true)(8).inputTransformer!((x, h) => OBFEval!(Model.BasisFunctions)(x, h))(coefs[6].dup);
    //auto st7c = new FIRState!(Complex!float, true)(8).inputTransformer!((x, h) => OBFEval!(Model.BasisFunctions)(x, h))(coefs[7].dup);
    

    //auto st5 = (new FIRState!(Complex!float, 16, true)).inputTransformer!((x, h) => OBFEval!(Model.BasisFunctions)(x, h))(coefs[4].dup);
    //auto st7 = (new FIRState!(Complex!float, 16, true)).inputTransformer!((x, h) => OBFEval!(Model.BasisFunctions)(x, h))(coefs[5].dup);

    //writeln("return filter");

    immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

    return filterBuilder(
            //st12,
            //makeRLSAdapter(st12, 1 - 1E-4, 1E-7),
            //st12c,
            //makeRLSAdapter(st12c, 1 - 1E-4, 1E-7),
            //st32,
            //makeRLSAdapter(st32, 1 - 1E-4, 1E-7),
            //st32c,
            //makeRLSAdapter(st32c, 1 - 1E-4, 1E-7),
            //st52,
            //makeRLSAdapter(st52, 1 - 1E-4, 1E-7),
            //st52c,
            //makeRLSAdapter(st52c, 1 - 1E-4, 1E-7),
            //st12,
            //makeRLSAdapter(st12, 0.98, 1E-7),
            //st32,
            //makeRLSAdapter(st32, 0.98, 1E-7),
            //st1c_2,
            //makeRLSAdapter(st1c_2, 0.98, 1E-7),
            //st3_2,
            //makeRLSAdapter(st3_2, 0.92, 1E-7),
            //st3c_2,
            //makeRLSAdapter(st3c_2, 0.999, 1E-7),
            st1,
            lsAdapter(st1, samplesOfOnePeriod)
                .trainingLimit(samplesOfOnePeriod * model.learningCount)
                .ignoreHeadSamples(samplesOfOnePeriod * model.learningCount * 0),
            //lmsAdapter(st1, 0.001, 1024, 0.5),
            //makeRLSAdapter(st1, 0.999, 1E2),
            st1c,
            lsAdapter(st1c, samplesOfOnePeriod)
                .trainingLimit(samplesOfOnePeriod * model.learningCount)
                .ignoreHeadSamples(samplesOfOnePeriod * model.learningCount * 1),
            //lmsAdapter(st1c, 0.001, 1024, 0.5),
            //makeRLSAdapter(st1c, 0.999, 1E2),
            //st2,
            //lsAdapter(st2, samplesOfOnePeriod),
            //lmsAdapter(st2, 0.002, 1024, 0.5),
            //makeRLSAdapter(st2, /*0.999*/1E2E-7),
            st3,
            lsAdapter(st3, samplesOfOnePeriod)
                .trainingLimit(samplesOfOnePeriod * model.learningCount)
                .ignoreHeadSamples(samplesOfOnePeriod * model.learningCount * 2),
            //lmsAdapter(st3, 0.001, 1024, 0.5),
            //makeRLSAdapter(st3, 0.999, 1E2),
            st3c,
            lsAdapter(st3c, samplesOfOnePeriod)
                .trainingLimit(samplesOfOnePeriod * model.learningCount)
                .ignoreHeadSamples(samplesOfOnePeriod * model.learningCount * 3),
            //lmsAdapter(st3c, 0.001, 1024, 0.5),
            //makeRLSAdapter(st3c, 0.999, 1E2),
            //st4,
            //lsAdapter(st4, samplesOfOnePeriod),
            //lmsAdapter(st4, 0.002, 1024, 0.5),
            //makeRLSAdapter(st4, /*0.999*/1E2E-7),
            //st4a,
            //lmsAdapter(st4a, 0.0005, 1024, 0.5),
            st5,
            lsAdapter(st5, samplesOfOnePeriod)
                .trainingLimit(samplesOfOnePeriod * model.learningCount)
                .ignoreHeadSamples(samplesOfOnePeriod * model.learningCount * 4),
            //lmsAdapter(st5, 0.001, 1024, 0.5),
            //makeRLSAdapter(st5, 0.999, 1E2),
            st5c,
            lsAdapter(st5c, samplesOfOnePeriod)
                .trainingLimit(samplesOfOnePeriod * model.learningCount)
                .ignoreHeadSamples(samplesOfOnePeriod * model.learningCount * 5),
            //lmsAdapter(st1, 0.001, 1024, 0.5),
            //makeRLSAdapter(st1, 1 - 1E-2, 1E-7),
            //lmsAdapter(st5c, 0.001, 1024, 0.5),
            //makeRLSAdapter(st5c, 0.999, 1E2),
            //st7,
            //lsAdapter(st7, 2000),
            //*********lmsAdapter(st7, 0.001, 1024, 0.5),
            //makeRLSAdapter(st7, 1 - 1E-6, 1E-7),
            //st12,
            //lsAdapter(st1, 2000),
            //lmsAdapter(st12, 0.010, 1024, 0.5),
            //makeRLSAdapter(st7, 0.9997, 1E-7),
            //*********st7c,
            //lsAdapter(st7c, 1000),
            //*********lmsAdapter(st7c, 0.001, 1024, 0.5),
            //makeRLSAdapter(st7c, 1 - 1E-6, 1E-7),
            //st5,
            //lmsAdapter(st5, 0.0085, 1024, 0.5),
            //makeRLSAdapter(st5, 1, 1E-7),
            //st3,
            //lmsAdapter(st3, 0.01, 1024, 0.5),
            //makeRLSAdapter(st3, 1, 1E-7),
            //st1,
            //lmsAdapter(st1, 0.01, 1024, 0.5),
            //makeRLSAdapter(st1, 1, 1E-7),
            );
}


auto makeParallelHammersteinWithDCMethodFilter(bool isOrthogonalized, Mod)(Mod mod, Model model)
{
    import dffdd.filter.polynomial;
    return makeCascadeHammersteinFilter!(isOrthogonalized, parallelFilter)(mod, model);
}


auto makeFrequencyHammersteinFilter(Model model)
{
    return new FrequencyHammersteinFilter!((i, s) => lsAdapter(s, model.learningSymbols).trainingLimit(model.learningSymbols * model.learningCount),
                                            BasisFunctions)(model.ofdm.subCarrierMap,
                                                            64,
                                                            model.ofdm.numOfFFT,
                                                            model.ofdm.numOfCP,
                                                            model.ofdm.scaleOfUpSampling);
}

