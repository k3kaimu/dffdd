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
import std.exception : enforce;

import dffdd.blockdiagram.adder;
import dffdd.blockdiagram.amplifier;
import dffdd.blockdiagram.decimator;
import dffdd.blockdiagram.filter;
import dffdd.blockdiagram.iqmixer;
import dffdd.blockdiagram.mod.ofdm;
import dffdd.blockdiagram.noise;
import dffdd.blockdiagram.quantizer;
import dffdd.blockdiagram.utils;
import dffdd.filter.filter;
import dffdd.filter.lms;
import dffdd.filter.ls;
import dffdd.filter.orthfreqpolyfil;
import dffdd.filter.orthogonalize;
import dffdd.filter.ph_dcm;
import dffdd.filter.primitives;
import dffdd.filter.rls;
import dffdd.filter.state;
import dffdd.filter.taylor;
import dffdd.utils.fft;
import dffdd.utils.unit;
import dffdd.utils.distribution;

import dffdd.mod.primitives;
import dffdd.mod.bpsk;
import dffdd.mod.qpsk;
import dffdd.mod.qam;
import dffdd.mod.ofdm;

import dffdd.dsp.statistics;


import constant;
import snippet;

enum size_t defaultDistortionOrder = 7;
alias CompleteDistorter(size_t P = defaultDistortionOrder) = PADistorter!(Complex!float, P);


struct Model
{
    size_t numOfModelTrainingSymbols = 300;
    size_t numOfFilterTrainingSymbols = 200;
    //size_t blockSize = 1024;
    size_t blockSize() const @property { return ofdm.numOfSamplesOf1Symbol*4; }
    real carrFreq = 2.45e9;
    real samplingFreq = 20e6 * 8;
    Gain SNR = 20.dB;
    Gain INR = 60.dB;
    bool withSIC = true;
    bool withSI = true;
    size_t learningSymbols = 10;
    size_t learningCount = 10;
    size_t swappedSymbols = 0;
    bool outputWaveform = false;
    bool outputBER = true;
    bool outputEVM = true;

    uint rndSeed = 114514;


    bool useDTXIQ = true;
    bool useDTXPN = false;
    bool useDTXPA = true;
    bool useSTXIQ = true;
    bool useSTXPN = false;
    bool useSTXPA = true;
    bool useSRXLN = true;
    bool useSRXIQ = true;
    bool useSRXQZ = true;

    bool doNoiseElimination = false;

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
        uint scaleOfUpSampling = 8;
        Gain PAPR = 10.dB;                 // 10dB
        //uint model.numOfSamplesOf1Symbol = 
        uint numOfSamplesOf1Symbol() const @property { return scaleOfUpSampling * (numOfFFT + numOfCP); }
        bool[] subCarrierMap() const @property
        {
            bool[] ret = new bool[scaleOfUpSampling * numOfFFT];

            // assert(numOfSubcarrier % 2 == 0);

            if(numOfSubcarrier % 2 == 0){
                ret[1 .. numOfSubcarrier/2 + 1] = true;
                ret[$ - numOfSubcarrier/2 .. $] = true;
            }else{
                ret[1 .. (numOfSubcarrier+1)/2 + 1] = true;
                ret[$ - (numOfSubcarrier+1)/2 + 1 .. $] = true;
            }

            return ret;
        }
    }
    OFDM ofdm;


    struct ThermalNoise
    {
        real temperature = 300;
        Voltage power(Model model)
        {
            return Voltage(sqrt(noisePower(model.samplingFreq, model.thermalNoise.temperature)));
        }
    }
    ThermalNoise thermalNoise;


    struct PowerSpectralDensityAnalyzer
    {
        uint resolution = 1024;
    }
    PowerSpectralDensityAnalyzer powerSpectralDensityAnalyzer;


    struct BERCounter
    {
        ulong totalBits = 10_000_000;
        ulong evmSymbols = 1000;
    }
    BERCounter berCounter;


    struct TXIQMixer
    {
        /* f(x) = x + b*x としたときの b */
        Complex!real imbCoef;
    }
    TXIQMixer txIQMixer;


    struct RXIQMixer
    {
        Complex!real imbCoef;
    }
    RXIQMixer rxIQMixer;


    struct PhaseNoise
    {
        real paramC = 1.65e-19;
    }
    PhaseNoise phaseNoise;


    struct PA 
    {
        // // https://datasheets.maximintegrated.com/en/ds/MAX2612-MAX2616.pdf
        // // MAX2616
        // Gain GAIN = 18.dB;
        // Voltage Vsat = (IIP3.dBm - 36).dB.gain.Voltage;
        // Voltage IIP3 = (37.2 - 18.4).dBm;
        // Voltage TX_POWER = 15.dBm;
        // https://datasheets.maximintegrated.com/en/ds/MAX2612-MAX2616.pdf
        // MAX2616
        Gain GAIN = 28.5.dB;
        Voltage Vsat = Voltage(21.8.dBm.volt / 2);
        Voltage IIP3 = 21.8.dBm;
        Voltage TX_POWER = 15.dBm;
        Gain MAX_VAR_IIP3 = 0.dB;       // IIP3の最大変位，dB単位で一様分布
        Gain MAX_VAR_TXP = 0.dB;        // 送信電力の最大変異，dB単位で一様分布
        Gain MAX_VAR_GAIN = 0.dB;       // 利得の最大変異，dB単位で一様分布
    }
    PA pa;

    struct LNA
    {
        Gain GAIN = 20.dB;      // 20dB
        Gain NF = 4.dB;         // 4dB
        uint noiseSeedOffset = 123;
        Gain DR = 70.dB;        // Dynamic Range
        Voltage IIP3 = (-3).dBm;  // MAX2695 https://datasheets.maximintegrated.com/en/ds/MAX2692-MAX2695.pdf
        uint smoothFactor = 3;
    }
    LNA lna;


    struct Quantizer
    {
        uint numOfBits = 14;
    }
    Quantizer quantizer;


    struct AWGN
    {
        // Gain SNR = 24;      // 24dB
    }
    AWGN awgn;


    //struct SIChannel
    //{
    //    size_t taps = 64;
    //    Gain c = (40.0 / 64).dB;     // 64サンプル遅延で80dB減衰
    //    bool isCoaxialCable = false;
    //}
    //SIChannel channel;
    struct SIChannel
    {
        Complex!real[] impulseResponse;
        size_t taps;
    }
    SIChannel channel;


    struct FIRFilter
    {
        size_t taps = 64;
    }
    FIRFilter firFilter;


    struct NLMSAdapter 
    {
        real mu;
    }
    NLMSAdapter nlmsAdapter;


    struct RLSAdapter
    {
        real delta;
        real lambda;
    }
    RLSAdapter rlsAdapter;


    struct LSAdapter {}
    LSAdapter lsAdapter;


    struct Orthogonalizer
    {
        bool enabled = false;
        size_t numOfTrainingSymbols = 10000;
    }
    Orthogonalizer orthogonalizer;


    struct BasisFunctionSelection
    {
        Gain noiseMargin = 6.dB;
        size_t nEstH = 2;
    }
    BasisFunctionSelection basisFuncsSelection;


    struct IterativeFreqSIC
    {
        size_t iterations = 2;
        size_t newtonIterations = 2;
    }
    IterativeFreqSIC iterativeFreqSIC;

    void txPower(Voltage p) @property
    {
        this.pa.TX_POWER = p;
    }


    void useCoaxialCableAsChannel() @property
    {
        this.channel.taps = 1;
    }
}


auto randomBits(uint prn, Model model)
{
    import dffdd.utils.binary;

    Random rnd;
    rnd.seed(prn * model.rndSeed);
    return dffdd.utils.binary.randomBits(rnd);
}


auto siRandomBits(Model model) { return randomBits(1, model); }
auto desiredRandomBits(Model model) { return randomBits(193, model); }



auto modOFDM(Model model)
{
    return chainedMod(
        //dffdd.mod.bpsk.BPSK.init,
        dffdd.mod.qam.QAM!(Complex!float)(16),
        new dffdd.mod.ofdm.OFDM!(Complex!float)(model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling),
    );
}


auto connectToModulator(R, Mod)(R r, Mod modObj, const(bool)* swExchange, Model model)
{
    auto swModObj = modObj.makeOFDMSymbolExchanger(swExchange, model.ofdm.numOfFFT, model.ofdm.numOfCP, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling);
    return new ModulatedRange!(R, typeof(swModObj))(r, swModObj);
}


auto connectToDemodulator(R, Mod)(R r, Mod modObj/*, Bits bits*/, Model)
{
    return new ModulatedRange!(R, Mod, true)(r, modObj);
}


auto thermalNoise(Model model, uint seedOffset = 123321)
{
    return ThermalNoise(model.samplingFreq, model.thermalNoise.temperature, model.rndSeed + seedOffset);
}


auto makeParallelHammersteinFilter(string optimizer, size_t distortionOrder = defaultDistortionOrder, bool useWL = true, bool isOrthogonalized, Mod)(Mod mod, Model model)
{
    alias C = Complex!float;

  static if(useWL)
    alias Dist = CompleteDistorter!(distortionOrder);
  else
  {
    // static assert(distortionOrder == 1);
    // alias Dist = Distorter!(C, x => x);
    alias Dist = OnlyPADistorter!(C, distortionOrder);
  }

    static if(isOrthogonalized)
    {
        alias GS = GramSchmidtOBFFactory!C;
        auto dist = new OrthogonalizedVectorDistorter!(C, Dist, GS)(new Dist(), new GS(Dist.outputDim));
    }
    else
    {
        auto dist = new Dist();
    }

    auto makeOptimizer(State)(State state)
    {
        immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

        static if(optimizer == "LMS")
            return makeNLMSAdapter(state, model.nlmsAdapter.mu).trainingLimit(samplesOfOnePeriod);
        else static if(optimizer == "RLS")
            return makeRLSAdapter(state, model.rlsAdapter.lambda, model.rlsAdapter.delta).trainingLimit(samplesOfOnePeriod);
        else static if(optimizer == "LS")
        {
            // immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;
            return makeLSAdapter(state, samplesOfOnePeriod).trainingLimit(samplesOfOnePeriod).ignoreHeadSamples(samplesOfOnePeriod);
        }
    }

    return new SimpleTimeDomainParallelHammersteinFilter!(Complex!float, typeof(dist), (s) => makeOptimizer(s))(dist, model.firFilter.taps);
}


auto makeCascadeHammersteinFilter(string optimizer, size_t distortionOrder = defaultDistortionOrder, Mod)(Mod mod, Model model)
{
    alias C = Complex!float;
    alias Dist = CompleteDistorter!(distortionOrder);
    alias GS = GramSchmidtOBFFactory!C;
    auto dist = new OrthogonalizedVectorDistorter!(C, Dist, GS)(new Dist(), new GS(Dist.outputDim));

    auto makeOptimizer(State)(size_t i, State state)
    {
        immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

        static if(optimizer == "LMS")
            return makeLMSAdapter(state, 0.02).trainingLimit(samplesOfOnePeriod);
        else static if(optimizer == "RLS")
            return makeRLSAdapter(state, 1, 1E-7).trainingLimit(samplesOfOnePeriod);
        else static if(optimizer == "LS")
        {
            return makeLSAdapter(state, 80 * 4 * model.learningSymbols).trainingLimit(samplesOfOnePeriod).ignoreHeadSamples(samplesOfOnePeriod * i);
        }
    }

    return new SimpleTimeDomainCascadeHammersteinFilter!(Complex!float, typeof(dist), (i, s) => makeOptimizer(i, s))(dist, model.firFilter.taps);
}


auto makeFrequencyHammersteinFilter2(string optimizer, size_t distortionOrder = defaultDistortionOrder, bool useWL = true)(Model model)
{
    import dffdd.filter.freqdomain;

    alias C = Complex!float;

    static if(useWL)
        alias Dist = CompleteDistorter!(distortionOrder);
    else
        alias Dist = OnlyPADistorter!(C, distortionOrder);

    auto makeOptimizer(State)(State state)
    {
        immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

      static if(optimizer == "LMS")
        return makeNLMSAdapter(state, model.nlmsAdapter.mu).trainingLimit(model.learningSymbols);
      else static if(optimizer == "RLS")
        return makeRLSAdapter(state, model.rlsAdapter.lambda, model.rlsAdapter.delta).trainingLimit(model.learningSymbols);
      else static if(optimizer == "LS")
        return makeLSAdapter(state, model.learningSymbols).trainingLimit(model.learningSymbols);
    }

    alias Adapter = typeof(makeOptimizer!(MultiFIRState!C)(MultiFIRState!C.init));

    auto dist = new Dist();
    auto regen = new OverlapSaveRegenerator2!C(Dist.outputDim, model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling);
    auto stateAdapter = new FrequencyDomainParallelHammersteinStateAdapter!(C, Adapter)
                        (
                            (size_t i, bool b, MultiFIRState!C s) => makeOptimizer(s),
                            model.ofdm.subCarrierMap,
                            Dist.outputDim
                        );

    return new FrequencyDomainHammersteinFilter!(
            Complex!float,
            typeof(dist),
            typeof(stateAdapter),
        )(
            dist,
            stateAdapter,
            model.ofdm.subCarrierMap,
            model.ofdm.numOfFFT,
            model.ofdm.numOfCP,
            model.ofdm.scaleOfUpSampling,
            model.samplingFreq,
            model.doNoiseElimination
        );
}


auto makeFrequencyDomainBasisFunctionSelector(Canceller)(Model model, Canceller canceller)
{
    import dffdd.filter.freqdomain;

    return new FrequencyDomainBasisFunctionSelector!(Complex!float, Canceller)(
        canceller,
        model.ofdm.subCarrierMap,
        model.ofdm.numOfFFT,
        model.ofdm.numOfCP,
        model.ofdm.scaleOfUpSampling,
        model.samplingFreq,
        model.basisFuncsSelection.nEstH,
        model.basisFuncsSelection.noiseMargin,
    );
}


auto makeFrequencyDCMHammersteinFilter(size_t type, Flag!"isParallel" isParallel, string optimizer, size_t distortionOrder = defaultDistortionOrder)(Model model)
if(type == 1 || type == 2)
{
    // alias BFs = BasisFunctions[0 .. numOfBasisFuncs];
    alias Dist = CompleteDistorter!(distortionOrder);

    auto makeOptimizer(State)(size_t p, State state)
    {
        immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

      static if(optimizer == "LMS")
        return makeLMSAdapter(state, 0.30).trainingLimit(model.learningSymbols);
      else static if(optimizer == "RLS")
        return makeRLSAdapter(state, 0.5, 1E-7).trainingLimit(model.learningSymbols);
      else static if(optimizer == "LS")
        return makeLSAdapter(state, model.learningSymbols).trainingLimit(model.learningCount * model.learningSymbols).ignoreHeadSamples(p * model.learningSymbols);
    }

    alias C = Complex!float;
    // alias Dist = Distorter!(C, BFs);

    auto dist = new Dist();

  static if(type == 1)
    alias FilterType = SimpleFrequencyDomainDCMHammersteinFilterType1;
  else
    alias FilterType = SimpleFrequencyDomainDCMHammersteinFilterType2; 

    return new FilterType!(
            Complex!float,
            typeof(dist),
            (i, bIsSC, p, s) => makeOptimizer(p, s),
            isParallel
        )(
            dist,
            model.ofdm.subCarrierMap,
            model.ofdm.numOfFFT,
            model.ofdm.numOfCP,
            model.ofdm.scaleOfUpSampling
        );
}


auto makeFrequencyDCMHammersteinFilter2(size_t type, Flag!"isParallel" isParallel, string optimizer, size_t distortionOrder = defaultDistortionOrder)(Model model, bool selectBF = false)
if(type == 2)
{
    import dffdd.filter.freqdomain;
    // alias BFs = BasisFunctions[0 .. numOfBasisFuncs];
    alias Dist = CompleteDistorter!(distortionOrder);

    auto makeOptimizer(State)(size_t p, State state)
    {
        immutable samplesOfOnePeriod = model.ofdm.numOfSamplesOf1Symbol * model.learningSymbols;

      static if(optimizer == "LMS")
        return makeLMSAdapter(state, 0.30).trainingLimit(model.learningSymbols);
      else static if(optimizer == "RLS")
        return makeRLSAdapter(state, 0.5, 1E-7).trainingLimit(model.learningSymbols);
      else static if(optimizer == "LS")
        return makeLSAdapter(state, model.learningSymbols).trainingLimit(model.learningCount * model.learningSymbols).ignoreHeadSamples(p * model.learningSymbols);
    }

    alias C = Complex!float;
    alias Adapter = typeof(makeOptimizer!(MultiFIRState!C)(0, MultiFIRState!C.init));

    auto dist = new Dist();
    auto regen = new OverlapSaveRegenerator2!C(Dist.outputDim, model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling);
    auto stateAdapter = new FrequencyDomainDCMHammersteinStateAdapter!(C, Adapter, isParallel)
                        (
                            (size_t i, bool b, size_t p, MultiFIRState!C s) => makeOptimizer(p, s),
                            model.ofdm.subCarrierMap,
                            Dist.outputDim
                        );

    return new FrequencyDomainHammersteinFilter!(
            Complex!float,
            typeof(dist),
            typeof(stateAdapter),
        )(
            dist,
            stateAdapter,
            model.ofdm.subCarrierMap,
            model.ofdm.numOfFFT,
            model.ofdm.numOfCP,
            model.ofdm.scaleOfUpSampling,
            model.samplingFreq,
            selectBF,
            model.basisFuncsSelection.nEstH,
            (-model.basisFuncsSelection.imageMargin.dB).dB,
            model.basisFuncsSelection.noiseMargin
        );
}


auto makeTaylorApproximationFilter(size_t distortionOrder = defaultDistortionOrder, size_t useWL = true)(Model model)
{
    alias C = Complex!float;

  static if(useWL)
    alias Dist = CompleteDistorter!(distortionOrder);
  else
  {
    static assert(distortionOrder == 1);
    alias Dist = Distorter!(C, x => x);
  }

    auto dist = new Dist();

    return new TaylorApproximationSICanceller!(Complex!float, typeof(dist))(dist, model.learningSymbols * model.ofdm.numOfSamplesOf1Symbol);
}
