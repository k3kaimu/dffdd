module models;

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

import dffdd.itpp;
import dffdd.mod.primitives;
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
                              //x => x * (x.re^^2 + x.im^^2)^^4,
                              //x => x.conj * (x.re^^2 + x.im^^2)^^4,
                              );

struct Model
{
    size_t numOfModelTrainingSample = 1024*1024;
    size_t numOfFilterTrainingSample = 1024*1024;
    size_t blockSize = 1024;
    real samplingFreq = 20e6 * 4;
    real SNR = 20;
    real INR = 60;
    bool withSIC = true;
    bool withSI = true;
    size_t numOfPopFront = 1;


    struct QAM
    {
        uint arity = 16;
    }
    QAM qam;


    struct OFDM
    {
        uint numOfFFT = 64;
        uint numOfCP = 16;
        uint numOfSubcarrier = 52;
        uint scaleOfUpSampling = 4;
        real PAPR = 10;                 // 10dB
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
        uint numOfBits = 14;
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
        dffdd.mod.qam.QAM(model.qam.arity),
        dffdd.mod.ofdm.OFDM(model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling),
    );
}


auto connectToModulator(R, Mod)(R r, Mod modObj, Model)
{
    return r
    .splitN(modObj.symInputLength)
    .tmap!(reverseArgs!mod, [0])(modObj)
    .joiner()
    .toWrappedRange
    ;
}


auto connectToDemodulator(R, Mod)(R r, Mod modObj/*, Bits bits*/, Model)
{
    return r
    .map!"cast(cfloat)a"
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


auto connectToTxChain(R)(R r)
{
    return r.connectToTXIQMixer().connectToPowerAmplifier();
}


auto connectToRxChain(R)(R r)
{
    return r.connectToLNA().connectToRXIQMixer().connectToQuantizer();
}


auto makeParallelHammersteinFilter(Mod)(Mod mod, Model model)
{
    import dffdd.filter.diagonal;
    import dffdd.filter.lms;
    import dffdd.filter.ls;
    import dffdd.filter.rls;
    import dffdd.filter.mempoly;
    import dffdd.filter.polynomial;
    import dffdd.filter.orthogonalize;

    auto state = new MemoryPolynomialState!(Complex!float, 8, 1, 0, 0, false, false)(1);

    //writeln("Use");
    //auto adapter = new LMSAdapter!(typeof(state))(state, 0.001, 1024, 0.5);
    //auto adapter = makeRLSAdapter(state, 1 - 1E-4, 1E-7);
    auto adapter = lsAdapter(state, 10000);

    return new PolynomialFilter!(typeof(state), typeof(adapter))(state, adapter);
}


auto makeCascadeHammersteinFilter(Mod)(Mod mod, Model model)
{
    import dffdd.filter.diagonal;
    import dffdd.filter.lms;
    import dffdd.filter.ls;
    import dffdd.filter.rls;
    import dffdd.filter.mempoly;
    import dffdd.filter.polynomial;
    import dffdd.filter.orthogonalize;
    import std.meta;
    import std.stdio;

    //writeln("orthogonalizer start");

    //auto orthogonalizer = new GramSchmidtOBFFactory!(Complex!float, BasisFunctions)();
    auto orthogonalizer = new DiagonalizationOBFFactory!(Complex!float, BasisFunctions)();

    {
        orthogonalizer.start();
        scope(exit)
            orthogonalizer.finish();


        auto signal = randomBits(1, model).connectToModulator(mod, model);

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

    Complex!float[][BasisFunctions.length] coefs;
    foreach(i, ref e; coefs){
        e = new Complex!float[BasisFunctions.length];
        orthogonalizer.getCoefs(i, e);
    }

    auto st1 = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[0].dup);
    auto st12 = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[0].dup);
    auto st1c = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[1].dup);
    auto st12c = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[1].dup);
    //auto st2 = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[2].dup);
    auto st3 = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[2].dup);
    auto st32 = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[2].dup);
    //auto st32 = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[2].dup);
    auto st3c = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[3].dup);
    auto st32c = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[3].dup);
    auto st5 = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[4].dup);
    auto st52 = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[4].dup);
    auto st5c = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[5].dup);
    auto st52c = (new FIRFilter!(Complex!float, 2, true)).inputTransformer!((x, h) => OBFEval!(BasisFunctions)(x, h))(coefs[5].dup);
    //auto st7 = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!(Model.BasisFunctions)(x, h))(coefs[6].dup);
    //auto st7c = (new FIRFilter!(Complex!float, 8, true)).inputTransformer!((x, h) => OBFEval!(Model.BasisFunctions)(x, h))(coefs[7].dup);
    

    //auto st5 = (new FIRFilter!(Complex!float, 16, true)).inputTransformer!((x, h) => OBFEval!(Model.BasisFunctions)(x, h))(coefs[4].dup);
    //auto st7 = (new FIRFilter!(Complex!float, 16, true)).inputTransformer!((x, h) => OBFEval!(Model.BasisFunctions)(x, h))(coefs[5].dup);

    //writeln("return filter");

    return serialFilter(
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
            lsAdapter(st1, 10000),
            //lmsAdapter(st1, 0.001, 1024, 0.5),
            //makeRLSAdapter(st1, 1 - 1E-4, 1E-7),
            st1c,
            lsAdapter(st1c, 10000),
            //lmsAdapter(st1c, 0.001, 1024, 0.5),
            //makeRLSAdapter(st1c, 1 - 1E-4, 1E-7),
            //st2,
            //lsAdapter(st2, 10000),
            //lmsAdapter(st2, 0.002, 1024, 0.5),
            //makeRLSAdapter(st2, /*0.9997*/1, 1E-7),
            st3,
            lsAdapter(st3, 10000),
            //lmsAdapter(st3, 0.001, 1024, 0.5),
            //makeRLSAdapter(st3, 1 - 1E-4, 1E-7),
            st3c,
            lsAdapter(st3c, 10000),
            //lmsAdapter(st3c, 0.001, 1024, 0.5),
            //makeRLSAdapter(st3c, 1 - 1E-4, 1E-7),
            //st4,
            //lsAdapter(st4, 10000),
            //lmsAdapter(st4, 0.002, 1024, 0.5),
            //makeRLSAdapter(st4, /*0.9997*/1, 1E-7),
            //st4a,
            //lmsAdapter(st4a, 0.0005, 1024, 0.5),
            st5,
            lsAdapter(st5, 10000),
            //lmsAdapter(st5, 0.001, 1024, 0.5),
            //makeRLSAdapter(st5, 1 - 1E-4, 1E-7),
            st5c,
            lsAdapter(st5c, 10000),
            //lmsAdapter(st5c, 0.001, 1024, 0.5),
            //makeRLSAdapter(st5c, 1 - 1E-4, 1E-7),
            //st7,
            //lsAdapter(st7, 10000),
            //*********lmsAdapter(st7, 0.001, 1024, 0.5),
            //makeRLSAdapter(st7, 1 - 1E-6, 1E-7),
            //st12,
            //lsAdapter(st1, 10000),
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