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
import dffdd.filter.lms;
import dffdd.filter.ls;
import dffdd.filter.orthfreqpolyfil;
import dffdd.filter.orthogonalize;
import dffdd.filter.ph_dcm;
import dffdd.filter.filter;
import dffdd.filter.primitives;
import dffdd.filter.rls;
import dffdd.filter.state;
import dffdd.utils.fft;
import dffdd.utils.unit;
//import dffdd.utils.msgpackrpc;




//import carbon.range;
import dranges.range;
import dranges.algorithm;

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
    size_t numOfFilterTrainingSymbols = 1000;
    //size_t blockSize = 1024;
    size_t blockSize() const @property { return ofdm.numOfSamplesOf1Symbol*4; }
    real carrFreq = 2.45e9;
    real samplingFreq = 20e6 * 4;
    real SNR = 20;
    real INR = 60;
    bool withSIC = true;
    bool withSI = true;
    size_t learningSymbols = 10;
    size_t learningCount = 10;
    size_t swappedSymbols = 0;
    bool outputWaveform = false;
    bool outputBER = true;


    bool useDTXIQ = true;
    bool useDTXPN = false;
    bool useDTXPA = true;
    bool useDTXIQ2 = false;
    bool useSTXIQ = true;
    bool useSTXPN = false;
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
            return noisePower(model.samplingFreq, model.thermalNoise.temperature);
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


    struct PhaseNoise
    {
        real paramC = 1.65e-19;
    }
    PhaseNoise phaseNoise;


    struct PA 
    {
        real GAIN = 27;     // 30dB
        real IIP3 = 17;     // 20dBm
        real TX_POWER = 15; // 20dBm
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
        uint numOfBits = 14;
    }
    Quantizer quantizer;


    struct AWGN
    {
        real SNR = 24;      // 24dB
    }
    AWGN awgn;


    struct SIChannel
    {
        size_t taps = 64;
    }
    SIChannel channel;


    struct FIRFilter
    {
        size_t taps = 64;
    }
    FIRFilter firFilter;


    struct Orthogonalizer
    {
        bool enabled = false;
        size_t numOfTrainingSymbols = 10000;
    }
    Orthogonalizer orthogonalizer;
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


auto connectToTXIQPhaseNoise(R)(R r, Model model)
{
    return r.connectTo!(dffdd.blockdiagram.iqmixer.PhaseNoise)(model.carrFreq, model.samplingFreq, model.phaseNoise.paramC);
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

    BoxMuller!Random gGen = BoxMuller!Random(rGen);

    Complex!float[] coefs;

    foreach(i; 0 .. 64){
        auto db = -40.0L * i / 64;
        coefs ~= cast(Complex!float)(gGen.front * 10.0L^^(db/20));
        gGen.popFront();
    }

    channelCoefsOfSI = coefs.dup;
}


auto connectToMultiPathChannel(R)(R r, Model model)
{
    return r.connectTo!FIRFilter(channelCoefsOfSI[0 .. model.channel.taps]).toWrappedRange;
}


auto connectToLNA(R)(R r, Model model)
{
    return r
    .add(thermalNoise(model).tmap!"a*b"(sqrt(model.lna.NF.dB.gain^^2 - 1)))
    .tmap!"a*b"(model.lna.GAIN.dB.gain)
    .toWrappedRange;
}


auto connectToQuantizer(R)(R r, Model model)
{
    return r
    .connectTo!PowerControlAmplifier((30 - model.ofdm.PAPR + 4.76).dBm)
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


//auto makeParallelHammersteinFilter(bool isOrthogonalized, string optimizer, size_t numOfBasisFuncs = BasisFunctions.length, Mod)(Mod mod, Model model)
//{
//    alias BFs = BasisFunctions[0 .. numOfBasisFuncs];

//    import dffdd.filter.lms;
//    import dffdd.filter.ls;
//    import dffdd.filter.rls;
//    import dffdd.filter.polynomial;
//    import dffdd.filter.orthogonalize;
//    import dffdd.filter.state;


//    auto orthogonalizer = new GramSchmidtOBFFactory!(Complex!float, BFs)();

//    {
//        orthogonalizer.start();
//        scope(exit)
//            orthogonalizer.finish();

//        auto signal = randomBits(1, model).connectToModulator(mod, new bool(false), model);
//        Complex!float[] buf = new Complex!float[1024];

//        foreach(i; 0 .. 1024){
//            foreach(j; 0 .. 1024){
//                auto f = signal.front;
//                buf[j] = complex(f.re, f.im);
//                signal.popFront();
//            }

//            .put(orthogonalizer, buf);
//        }
//    }

//    Complex!float[][BFs.length] coefs;
//    foreach(i, ref e; coefs){
//        e = new Complex!float[BFs.length];
//        orthogonalizer.getCoefs(i, e);
//    }


//    Complex!float delegate(Complex!float x)[BFs.length] bflist;
//    foreach(i, BF; BFs){
//      static if(isOrthogonalized)
//        bflist[i] = delegate Complex!float (Complex!float x) { return OBFEval!BFs(x, coefs[i]); };
//      else
//        bflist[i] = delegate Complex!float (Complex!float x) { return BF(x); };
//    }

//    auto state = new ParallelHammersteinState!(Complex!float, BFs.length, true)(model.firFilter.taps, bflist);

//  static if(optimizer == "LMS")
//    auto adapter = new LMSAdapter!(typeof(state))(state, 0.02);
//  else static if(optimizer == "RLS")
//    auto adapter = makeRLSAdapter(state, 1 - 1E-4, 1E-7);
//  else static if(optimizer == "LS")
//  {
//    immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;
//    auto adapter = lsAdapter(state, 80 * 4 * model.learningSymbols).trainingLimit(samplesOfOnePeriod * model.learningCount).ignoreHeadSamples(samplesOfOnePeriod);
//  }

//    return oneStateFilter(state, adapter);


//    /+  GeneralParallelHammersteinFilterを使う場合はこんな感じになる(コメントアウト最後まで)
//    Complex!float[BFs.length][] distortionFunc(in Complex!float[] tx)
//    {
//        Complex!float[BFs.length][] dst;

//        foreach(e; tx){
//            Complex!float[BFs.length] ds;

//            foreach(i, BF; BFs){
//              static if(isOrthogonalized)
//                ds[i] = OBFEval!BFs(e, coefs[i]);
//              else
//                ds[i] = BF(x);
//            }

//            dst ~= ds;
//        }

//        return dst;
//    }

//  static if(optimizer == "LMS")
//    auto makeAdapter(State)(State state){ return lmsAdapter(state, 0.02); }
//  else static if(optimizer == "RLS")
//    auto makeAdapter(State)(State state){ return makeRLSAdapter(state, 1 - 1E-4, 1E-7); }
//  else static if(optimizer == "LS")
//  {
//    immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;
//    auto makeAdapter(State)(State state){ return lsAdapter(state, 80 * 4 * model.learningSymbols).trainingLimit(samplesOfOnePeriod * model.learningCount).ignoreHeadSamples(samplesOfOnePeriod); }
//  }

//    return generalParallelHammersteinFilter!(Complex!float, BFs.length, distortionFunc, makeAdapter)(model.firFilter.taps);
//    +/
//}


///**
//Fast Training and Cancelling method for Hammerstein Self-Interference Canceller on Full-Duplex OFDM Communication
//*/
//auto makeCascadeHammersteinFilterImpl(bool isOrthogonalized,
//                                  string optimizer,
//                                  alias filterBuilder = serialFilter,
//                                  bool isSerialized = true,
//                                  alias makeStructure = (n) => iota(n).map!"[a]",
//                                  Mod)
//                                (Mod mod, Model model)
//{
//    import dffdd.filter.lms;
//    import dffdd.filter.ls;
//    import dffdd.filter.rls;
//    import dffdd.filter.polynomial;
//    import dffdd.filter.orthogonalize;
//    import dffdd.filter.state;
//    import std.meta;
//    import std.stdio;

//    //writeln("orthogonalizer start");

//    auto orthogonalizer = new GramSchmidtOBFFactory!(Complex!float, BasisFunctions)();
//    //auto orthogonalizer = new DiagonalizationOBFFactory!(Complex!float, BasisFunctions)();

//    {
//        orthogonalizer.start();
//        scope(exit)
//            orthogonalizer.finish();

//        auto signal = randomBits(1, model).connectToModulator(mod, new bool(false), model);

//        //.put(orthogonalizer, signal.take(1024*400));
//        Complex!float[] buf = new Complex!float[1024];
//        foreach(i; 0 .. 1024){
//            foreach(j; 0 .. 1024){
//                auto f = signal.front;
//                buf[j] = complex(f.re, f.im);
//                signal.popFront();
//            }
//            .put(orthogonalizer, buf);
//        }
//    }


//  static if(isOrthogonalized)
//  {
//    pragma(msg, "Orthogonalized");
//    Complex!float[][BasisFunctions.length] coefs;
//    foreach(i, ref e; coefs){
//        e = new Complex!float[BasisFunctions.length];
//        orthogonalizer.getCoefs(i, e);
//    }
//  }
//  else
//  {
//    pragma(msg, "Un-Orthogonalized");
//  }

//    static
//    string makeDefFilter()
//    {
//        auto app = appender!string();
//        auto braches = makeStructure(BasisFunctions.length).map!"a.array()".array();

//        // 各ブランチの基底関数の宣言
//        foreach(i, blist; braches){
//            app.formattedWrite("Complex!float delegate(Complex!float x)[%2$s] bflist%1$s;", i, blist.length);
//            foreach(j, e; blist){
//                if(isOrthogonalized)
//                    app.formattedWrite("bflist%1$s[%2$s] = delegate Complex!float (Complex!float x) { return OBFEval!(BasisFunctions[%3$s])(x, coefs[%3$s]); };\n", i, j, e);
//                else
//                    app.formattedWrite("bflist%1$s[%2$s] = delegate Complex!float (Complex!float x) { return BasisFunctions[%3$s](x); };\n", i, j, e);
//            }
//        }

//        // 各ブランチの状態の宣言
//        foreach(i, blist; braches){
//            if(blist.length != 1)
//                app.formattedWrite("auto st%1$s = new ParallelHammersteinState!(Complex!float, %2$s, true)(model.firFilter.taps, bflist%1$s);\n", i, blist.length);
//            else
//                app.formattedWrite("auto st%1$s = new FIRState!(Complex!float, true)(model.firFilter.taps).inputTransformer!((x, f) => f(x))(bflist%1$s[0]);\n", i);
//        }

//        // 全部ランチ結合
//        app.formattedWrite("auto result = filterBuilder(");
//        foreach(i, blist; braches)
//            app.formattedWrite("st%1$s, makeOptimizer!(%1$s + 1)(st%1$s, model), ", i);
//        app.formattedWrite(");");
//        return app.data;
//    }

//    static
//    auto makeOptimizer(size_t p, State)(State state, const ref Model model)
//    {
//        immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

//      static if(optimizer == "LMS")
//        return lmsAdapter(state, 0.01);
//      else static if(optimizer == "RLS")
//        return makeRLSAdapter(state, 0.999, 1E-7);
//      else static if(optimizer == "LS")
//        return lsAdapter(state, samplesOfOnePeriod)
//                .trainingLimit(samplesOfOnePeriod * model.learningCount)
//                .ignoreHeadSamples(isSerialized ? samplesOfOnePeriod * model.learningCount * (p-1) + samplesOfOnePeriod * p : samplesOfOnePeriod);
//    }

//    mixin(makeDefFilter());
//    return result;
//}


//auto makeCascadeHammersteinFilter(bool isOrthogonalized, string optimizer, alias filterBuilder = serialFilter, bool isSerialized = true, Mod)(Mod mod, Model model)
//{
//    return makeCascadeHammersteinFilterImpl!(isOrthogonalized, optimizer, filterBuilder, isSerialized, (n) => iota(n).map!"[a]")(mod, model);
//}


//auto makeCascadeWLHammersteinFilter(bool isOrthogonalized, string optimizer, Mod)(Mod mod, Model model)
//{
//    return makeCascadeHammersteinFilterImpl!(isOrthogonalized, optimizer, serialFilter, true, (n) => iota(n).chunks(2))(mod, model);
//}


//auto makeCascadeWL1HammersteinFilter(bool isOrthogonalized, string optimizer, Mod)(Mod mod, Model model)
//{
//    return makeCascadeHammersteinFilterImpl!(isOrthogonalized, optimizer, serialFilter, true, (n) => [0, 1] ~ iota(2, n).map!"[a]".array())(mod, model);
//}


//auto makeParallelHammersteinWithDCMethodFilter(bool isOrthogonalized, string optimizer = "LS", Mod)(Mod mod, Model model)
//{
//    return makeCascadeHammersteinFilter!(isOrthogonalized, optimizer, parallelFilter, false)(mod, model);


//  //  auto orthogonalizer = new GramSchmidtOBFFactory!(Complex!float, BasisFunctions)();
//  //  //auto orthogonalizer = new DiagonalizationOBFFactory!(Complex!float, BasisFunctions)();

//  //  {
//  //      orthogonalizer.start();
//  //      scope(exit)
//  //          orthogonalizer.finish();

//  //      auto signal = randomBits(1, model).connectToModulator(mod, new bool(false), model);

//  //      Complex!float[] buf = new Complex!float[1024];
//  //      foreach(i; 0 .. 1024){
//  //          foreach(j; 0 .. 1024){
//  //              auto f = signal.front;
//  //              buf[j] = complex(f.re, f.im);
//  //              signal.popFront();
//  //          }
//  //          .put(orthogonalizer, buf);
//  //      }
//  //  }


//  //static if(isOrthogonalized)
//  //{
//  //  pragma(msg, "Orthogonalized");
//  //  Complex!float[][BasisFunctions.length] coefs;
//  //  foreach(i, ref e; coefs){
//  //      e = new Complex!float[BasisFunctions.length];
//  //      orthogonalizer.getCoefs(i, e);
//  //  }


//  //  auto adaptInputTransformer(size_t i, S)(S state)
//  //  {
//  //      return state.inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[i].dup);
//  //  }
//  //}
//  //else
//  //{
//  //  pragma(msg, "Un-Orthogonalized");
//  //  auto adaptInputTransformer(size_t i, S)(S state)
//  //  {
//  //      return state.inputTransformer!(x => BasisFunctions[i](x))();
//  //  }
//  //}

//  //  auto st1 = adaptInputTransformer!0(new FIRState!(Complex!float, true)(model.firFilter.taps));
//  //  auto st1c = adaptInputTransformer!1(new FIRState!(Complex!float, true)(model.firFilter.taps));
//  //  auto st3 = adaptInputTransformer!2(new FIRState!(Complex!float, true)(model.firFilter.taps));
//  //  auto st3c = adaptInputTransformer!3(new FIRState!(Complex!float, true)(model.firFilter.taps));
//  //  auto st5 = adaptInputTransformer!4(new FIRState!(Complex!float, true)(model.firFilter.taps));
//  //  auto st5c = adaptInputTransformer!5(new FIRState!(Complex!float, true)(model.firFilter.taps));


//  //  return makeParallelHammersteinWithDCM
//  //          (model.firFilter.taps, model.learningSymbols * model.ofdm.numOfSamplesOf1Symbol, model.learningSymbols * model.learningCount,
//  //              st1, st1c, st3, st3c, st5, st5c);
//}


auto makeParallelHammersteinFilter(string optimizer, size_t numOfBasisFuncs = BasisFunctions.length, Mod)(Mod mod, Model model)
{


    alias BFs = BasisFunctions[0 .. numOfBasisFuncs];

    alias C = Complex!float;
    alias Dist = Distorter!(C, BFs);
    alias GS = GramSchmidtOBFFactory!C;
    auto dist = new OrthogonalizedVectorDistorter!(C, Dist, GS)(new Dist(), new GS(numOfBasisFuncs));

    auto makeOptimizer(State)(State state)
    {
        immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

        static if(optimizer == "LMS")
            return makeLMSAdapter(state, 0.02).trainingLimit(samplesOfOnePeriod);
        else static if(optimizer == "RLS")
            return makeRLSAdapter(state, 1 - 1E-4, 1E-7).trainingLimit(samplesOfOnePeriod);
        else static if(optimizer == "LS")
        {
            immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;
            return makeLSAdapter(state, samplesOfOnePeriod).trainingLimit(samplesOfOnePeriod).ignoreHeadSamples(samplesOfOnePeriod);
        }
    }

    return new SimpleTimeDomainParallelHammersteinFilter!(Complex!float, typeof(dist), (s) => makeOptimizer(s))(dist, numOfBasisFuncs, model.firFilter.taps);
}


auto makeCascadeHammersteinFilter(string optimizer, size_t numOfBasisFuncs = BasisFunctions.length, Mod)(Mod mod, Model model)
{
    alias BFs = BasisFunctions[0 .. numOfBasisFuncs];

    alias C = Complex!float;
    alias Dist = Distorter!(C, BFs);
    alias GS = GramSchmidtOBFFactory!C;
    auto dist = new OrthogonalizedVectorDistorter!(C, Dist, GS)(new Dist(), new GS(numOfBasisFuncs));

    auto makeOptimizer(State)(size_t i, State state)
    {
        immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

        static if(optimizer == "LMS")
            return makeLMSAdapter(state, 0.02).trainingLimit(samplesOfOnePeriod);
        else static if(optimizer == "RLS")
            return makeRLSAdapter(state, 1 - 1E-4, 1E-7).trainingLimit(samplesOfOnePeriod);
        else static if(optimizer == "LS")
        {
            return makeLSAdapter(state, 80 * 4 * model.learningSymbols).trainingLimit(samplesOfOnePeriod).ignoreHeadSamples(samplesOfOnePeriod * i);
        }
    }

    return new SimpleTimeDomainCascadeHammersteinFilter!(Complex!float, typeof(dist), (i, s) => makeOptimizer(i, s))(dist, numOfBasisFuncs, model.firFilter.taps);
}


auto makeFrequencyHammersteinFilter(string optimizer, size_t numOfBasisFuncs = BasisFunctions.length)(Model model)
{
    alias BFs = BasisFunctions[0 .. numOfBasisFuncs];

   auto makeOptimizer(State)(State state)
   {
       immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

     static if(optimizer == "LMS")
       return makeLMSAdapter(state, 0.30).trainingLimit(model.learningSymbols);
     else static if(optimizer == "RLS")
       return makeRLSAdapter(state, 0.97, 1E-7).trainingLimit(model.learningSymbols);
     else static if(optimizer == "LS")
       return makeLSAdapter(state, model.learningSymbols).trainingLimit(model.learningSymbols);
   }

   //auto distorter = Distorter!BasisFunctions();
    // auto distorter = new Distorter!(Complex!float, BasisFunctions)();
   //auto specConv = new BasisSpectrumConverter!(Complex!float, BasisFunctions.length)(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling);

    alias C = Complex!float;
    alias Dist = Distorter!(C, BFs);
    alias GS = GramSchmidtOBFFactory!C;
    auto dist = new OrthogonalizedVectorDistorter!(C, Dist, GS)(new Dist(), new GS(numOfBasisFuncs));

    auto specConv = new OrthogonalizedSpectrumConverter!(C, GS)(new GS(numOfBasisFuncs), numOfBasisFuncs, model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling);

   return new SimpleFrequencyDomainParallelHammersteinFilter!(
           //(i, bIsSC, s) => lsAdapter(s, model.learningSymbols).trainingLimit(model.learningSymbols * model.learningCount),
           //(i, bIsSC, s) => makeRLSAdapter(s, 0.97, 1E-7),//.trainingLimit(model.learningSymbols * model.learningCount),
           //(i, bIsSC, s) => lmsAdapter(s, 0.4, 1024, 0.5),
           Complex!float,
           typeof(dist),
           typeof(specConv),
           (i, bIsSC, s) => makeOptimizer(s),
        //    typeof(specConv),
       )(
           dist,
           specConv,
           model.ofdm.subCarrierMap,
        //    model.firFilter.taps,
           model.ofdm.numOfFFT,
           model.ofdm.numOfCP,
           model.ofdm.scaleOfUpSampling
       );
}


auto makeFrequencyCascadeHammersteinFilter(string optimizer, size_t numOfBasisFuncs = BasisFunctions.length)(Model model)
{
    alias BFs = BasisFunctions[0 .. numOfBasisFuncs];

   auto makeOptimizer(State)(State state)
   {
       immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

     static if(optimizer == "LMS")
       return makeLMSAdapter(state, 0.30).trainingLimit(model.learningSymbols);
     else static if(optimizer == "RLS")
       return makeRLSAdapter(state, 0.5, 1E-7).trainingLimit(model.learningSymbols);
     else static if(optimizer == "LS")
       return makeLSAdapter(state, model.learningSymbols).trainingLimit(model.learningSymbols);
   }

   //auto distorter = Distorter!BasisFunctions();
    // auto distorter = new Distorter!(Complex!float, BasisFunctions)();
   //auto specConv = new BasisSpectrumConverter!(Complex!float, BasisFunctions.length)(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling);

    alias C = Complex!float;
    alias Dist = Distorter!(C, BFs);
    alias GS = GramSchmidtOBFFactory!C;
    auto dist = new OrthogonalizedVectorDistorter!(C, Dist, GS)(new Dist(), new GS(numOfBasisFuncs));

    auto specConv = new OrthogonalizedSpectrumConverter!(C, GS)(new GS(numOfBasisFuncs), numOfBasisFuncs, model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling);

   return new SimpleFrequencyDomainCascadeHammersteinFilter!(
           //(i, bIsSC, s) => lsAdapter(s, model.learningSymbols).trainingLimit(model.learningSymbols * model.learningCount),
           //(i, bIsSC, s) => makeRLSAdapter(s, 0.97, 1E-7),//.trainingLimit(model.learningSymbols * model.learningCount),
           //(i, bIsSC, s) => lmsAdapter(s, 0.4, 1024, 0.5),
           Complex!float,
           typeof(dist),
           typeof(specConv),
           (i, bIsSC, p, s) => makeOptimizer(s),
        //    typeof(specConv),
       )(
           dist,
           specConv,
           model.ofdm.subCarrierMap,
        //    model.firFilter.taps,
           model.ofdm.numOfFFT,
           model.ofdm.numOfCP,
           model.ofdm.scaleOfUpSampling
       );
}


/+
auto makeFrequencyHammersteinFilter(string optimizer)(Model model)
{
    auto makeOptimizer(State)(State state)
    {
      static if(optimizer == "LMS")
        return makeLMSAdapter(state, 0.30);
      else static if(optimizer == "RLS")
        return makeRLSAdapter(state, 0.97, 1E-7);
      else static if(optimizer == "LS")
        return makeLSAdapter(state, model.learningSymbols).trainingLimit(model.learningSymbols * model.learningCount);
    }

    return new FrequencyHammersteinFilter!(
            //(i, bIsSC, s) => lsAdapter(s, model.learningSymbols).trainingLimit(model.learningSymbols * model.learningCount),
            //(i, bIsSC, s) => makeRLSAdapter(s, 0.97, 1E-7),//.trainingLimit(model.learningSymbols * model.learningCount),
            //(i, bIsSC, s) => lmsAdapter(s, 0.4, 1024, 0.5),
            (i, bIsSC, s) => makeOptimizer(s),
            BasisFunctions
        )(
            model.ofdm.subCarrierMap,
            model.firFilter.taps,
            model.ofdm.numOfFFT,
            model.ofdm.numOfCP,
            model.ofdm.scaleOfUpSampling
        );
}
+/


//auto makeAliasRemovableParallelHammersteinFilter(bool isOrthogonalized, string optimizer = "LS", Mod)(Mod mod, Model model)
//{
//    static assert(!isOrthogonalized, "Orthogonalization of ARPH have not implemented yet.");

//    immutable numOfFFT = model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling;
//    immutable numOfCP = model.ofdm.numOfCP * model.ofdm.scaleOfUpSampling;

//    auto distortionFunc(in Complex!float[] tx)
//    {
//        return generateOFDMAliasSignal!(1, BasisFunctions)(tx, numOfFFT, numOfCP);
//    }

//    enum size_t numOfFIRFilters = typeof(distortionFunc(null)[0]).init.length;

//  static if(optimizer == "LMS")
//    auto makeAdapter(State)(State state){ return lmsAdapter(state, 0.02); }
//  else static if(optimizer == "RLS")
//    auto makeAdapter(State)(State state){ return makeRLSAdapter(state, 1 - 1E-4, 1E-7); }
//  else static if(optimizer == "LS")
//  {
//    immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;
//    auto makeAdapter(State)(State state){ return lsAdapter(state, 80 * 4 * model.learningSymbols).trainingLimit(samplesOfOnePeriod * model.learningCount).ignoreHeadSamples(samplesOfOnePeriod); }
//  }

//    return generalParallelHammersteinFilter!(Complex!float, numOfFIRFilters, distortionFunc, makeAdapter)(model.firFilter.taps);
//}


// auto makeParallelHammersteinFilter(bool isOrthogonalized, string optimizer, size_t numOfBasisFuncs = BasisFunctions.length, Mod)(Mod mod, Model model)
// {
//    alias BFs = BasisFunctions[0 .. numOfBasisFuncs];

//    import dffdd.filter.lms : makeLMSAdapter;
//    import dffdd.filter.ls : makeLSAdapter;
//    import dffdd.filter.rls : makeRLSAdapter;
//    import dffdd.filter.polynomial : ParallelHammersteinState;
//    import dffdd.filter.orthogonalize : GramSchmidtOBFFactory;
//    import dffdd.filter.state;


//    auto orthogonalizer = new GramSchmidtOBFFactory!(Complex!float, BFs)();

//    {
//        orthogonalizer.start();
//        scope(exit)
//            orthogonalizer.finish();

//        auto signal = randomBits(1, model).connectToModulator(mod, new bool(false), model);
//        Complex!float[] buf = new Complex!float[1024];

//        foreach(i; 0 .. 1024){
//            foreach(j; 0 .. 1024){
//                auto f = signal.front;
//                buf[j] = complex(f.re, f.im);
//                signal.popFront();
//            }

//            .put(orthogonalizer, buf);
//        }
//    }

//    Complex!float[][BFs.length] coefs;
//    foreach(i, ref e; coefs){
//        e = new Complex!float[BFs.length];
//        orthogonalizer.getCoefs(i, e);
//    }


//    Complex!float delegate(Complex!float x)[BFs.length] bflist;
//    foreach(i, BF; BFs){
//      static if(isOrthogonalized)
//        bflist[i] = delegate Complex!float (Complex!float x) { return OBFEval!BFs(x, coefs[i]); };
//      else
//        bflist[i] = delegate Complex!float (Complex!float x) { return BF(x); };
//    }

//    auto state = new ParallelHammersteinState!(Complex!float, BFs.length, true)(model.firFilter.taps, bflist);

//  static if(optimizer == "LMS")
//    auto adapter = new LMSAdapter!(typeof(state))(state, 0.02);
//  else static if(optimizer == "RLS")
//    auto adapter = makeRLSAdapter(state, 1 - 1E-4, 1E-7);
//  else static if(optimizer == "LS")
//  {
//    immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;
//    auto adapter = lsAdapter(state, 80 * 4 * model.learningSymbols).trainingLimit(samplesOfOnePeriod * model.learningCount).ignoreHeadSamples(samplesOfOnePeriod);
//  }

//    return oneStateFilter(state, adapter);
// }
