module app;

import std.algorithm;
import std.complex;
import std.conv : text;
import std.math;
import std.random;
import std.range;
import std.stdio;
import std.traits;

import dffdd.mod;
import dffdd.utils.fft;
import dffdd.mod.sefdm;

import tuthpc.taskqueue;

alias F = float;
alias C = Complex!F;


extern(C) void openblas_set_num_threads(int num_threads);
extern(C) int openblas_get_num_threads();

shared static this()
{
    openblas_set_num_threads(1);
}


void main()
{
    import std.file : mkdirRecurse;

    auto env = defaultJobEnvironment();
    auto taskList = new MultiTaskList!void();

    static void run(alias simImpl, ParamType)(string filename, ParamType[] params, Flag!"shortcut" shortcut = No.shortcut)
    {
        import std.format;
        import std.json;
        import std.file;
        JSONValue[] berList;

        if(std.file.exists(filename))
            return;

        foreach(p; params) {
            auto simResult = simImpl(p);
            writefln!"\t%s: %s, MOD: %s Mbps, DEMOD: %s Mbps"(p.vEbN0dB, simResult.ber, simResult.modThroughput / 1e6, simResult.demodThroughput / 1e6);

            berList ~= JSONValue(simResult.ber);

            if(shortcut && simResult.ber < 1e-7)
                break;
        }
        writeln();

        std.file.write(filename, JSONValue(berList).toString());
    }

    // // AWGN
    // foreach(nFFT; [64, 128, 256, 512, 1024, 2048]) {
    //     mkdirRecurse("AWGN");
    //     foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5]) {
    //         taskList.append(&run, "AWGN", nFFT, nTone, true, 1, 20, Yes.isPerfectChEst, 1, "Perfect");
    //     }
    // }


    // // Rayleigh
    // foreach(nFFT; [64, 128, 256, 512, 1024, 2048]) {
    //     mkdirRecurse("Rayleigh");
    //     foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5])
    //     foreach(nChTap; [1, 4, 8, 16, 32, 64]) {
    //         taskList.append(&run, "Rayleigh", nFFT, nTone, false, nChTap, 20, Yes.isPerfectChEst, 1, "Perfect");
    //     }
    // }

    // foreach(nFFT; [2048]) {
    //     mkdirRecurse("Rayleigh_Iter");
    //     foreach(nIter; iota(2, 41, 2)) {
    //         taskList.append(&run, "Rayleigh_Iter", nFFT, nFFT / 8 * 5, false, 8, nIter, Yes.isPerfectChEst, 1, "Perfect");
    //     }
    // }

    foreach(chEstMethod; [ChannelEstimationMethod.perfect, ChannelEstimationMethod.ZF, ChannelEstimationMethod.MMSE]) {
        auto numSymChEstList = (chEstMethod == ChannelEstimationMethod.perfect) ? [1] : [1, 2, 4, 8, 16];

        // AWGN
        foreach(numSymChEst; numSymChEstList) {
            foreach(nFFT; [64, 128, 256, 512, 1024, 2048, 4096]) {
                mkdirRecurse("AWGN");
                foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5]) {
                    alias M = QPSK!C;

                    SimParams!M[] params;
                    foreach(vEbN0dB; iota(21)) {
                        SimParams!M p;
                        p.sefdm.nFFT = nFFT;
                        p.sefdm.nTone = nTone;
                        p.sefdm.nSpreadIn = 4096;
                        p.sefdm.nSpreadOut = 4096 / nFFT * nTone;
                        // p.sefdm.nData = nFFT;
                        p.sefdm.nIter = 20;
                        p.sefdm.rndSeed = 0;
                        p.vEbN0dB = vEbN0dB;
                        p.channelType = "AWGN";
                        p.channel = makeAWGNChannel!C(sigma2: 10.0^^(-vEbN0dB / 10.0) / M.symInputLength, seed: 0);
                        p.numSymChEst = numSymChEst;
                        p.chEstMethod = chEstMethod;
                        p.rndSeed = 0;

                        params ~= p;
                    }

                    auto filename = i"AWGN/nFFT$(nFFT)_nTone$(nTone)_nIter20_EstCh$(numSymChEst)_$(chEstMethod).json".text();
                    taskList.append(&run!(runSimImpl!(C, M), SimParams!M), filename, params, Yes.shortcut);
                }
            }
        }

        // Rayleigh
        foreach(numSymChEst; numSymChEstList) {
            foreach(nFFT; [64, 128, 256, 512, 1024, 2048, 4096]) {
                mkdirRecurse("Rayleigh");
                foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5])
                foreach(nChTap; [1, 2, 4, 8, 16, 32, 64, 128])
                {
                    if(nChTap > nFFT / 4)
                        continue;

                    alias M = QPSK!C;
                    immutable nIter = 20;

                    SimParams!M[] params;
                    foreach(vEbN0dB; iota(21)) {
                        SimParams!M p;
                        p.sefdm.nFFT = nFFT;
                        p.sefdm.nTone = nTone;
                        p.sefdm.nSpreadIn = 4096;
                        p.sefdm.nSpreadOut = 4096 / nFFT * nTone;
                        // p.sefdm.nData = nFFT;
                        p.sefdm.nIter = nIter;
                        p.sefdm.rndSeed = 0;
                        p.vEbN0dB = vEbN0dB;
                        p.channelType = "Rayleigh";
                        p.channel = makeRayleighChannel!C(nTaps: nChTap / 8, sigma2: 10.0^^(-vEbN0dB / 10.0) / M.symInputLength, seed: 0);
                        p.numSymChEst = numSymChEst;
                        p.chEstMethod = chEstMethod;
                        p.rndSeed = 0;

                        params ~= p;
                    }

                    auto filename = i"Rayleigh/nFFT$(nFFT)_nTone$(nTone)_nChTap$(nChTap)_nIter$(nIter)_EstCh$(numSymChEst)_$(chEstMethod).json".text();
                    taskList.append(&run!(runSimImpl!(C, M), SimParams!M), filename, params, Yes.shortcut);
                }
            }
        }
    }


    // foreach(nFFT; [64, 128, 256, 512, 1024, 2048]) {
    //     foreach(nTone; [/*nFFT, nFFT / 8 * 7, nFFT / 8 * 6,*/ nFFT / 8 * 5])
    //     foreach(nChTap; [1, /*nFFT / 32, nFFT / 32 * 2, nFFT / 32 * 3,*/ nFFT / 32 * 4]) {
    //         taskList.append(&run, nFFT, nTone, false, nChTap, 100);
    //     }
    // }


    tuthpc.taskqueue.run(taskList, env);
    // foreach(i; 0 .. taskList.length)
    //     taskList[i]();
}

struct SimResult
{
    double ber;
    long errbits;
    long totalbits;
    double modThroughput, demodThroughput;
}


struct RDFTsSEFDMParams(Mod)
{
    Mod mod;
    uint nFFT;
    uint nTone;
    uint nSpreadIn;
    uint nSpreadOut;
    uint rndSeed = 0;
    uint nIter;
    // alias nData = nFFT;
}


enum ChannelEstimationMethod : string
{
    perfect = "Perfect",
    ZF = "ZF",
    MMSE = "MMSE",
}


struct SimParams(Mod)
{
    RDFTsSEFDMParams!Mod sefdm;
    string channelType;
    WirelessChannel!C channel;
    double vEbN0dB;
    size_t maxErrbits = 1000;
    size_t maxTotalbits = 100_000_000;
    size_t minTotalbits = 1_000_000;
    // size_t nChTaps = 8;
    // bool isChNorm = false;
    ChannelEstimationMethod chEstMethod = ChannelEstimationMethod.ZF;
    size_t numSymChEst = 1;
    size_t rndSeed;
}


SimResult runSimImpl(C, Mod)(SimParams!Mod params)
{
    import std.datetime.stopwatch;

    // immutable double SIGMA2 = 1.0 / (10^^(params.vEbN0dB/10.0) * params.sefdm.mod.symInputLength) * params.sefdm.nFFT / params.sefdm.nData;

    size_t totalbits = 0;
    size_t errbits = 0;

    size_t rndSeed = 0;

    // auto sefdm = new RDFTsSEFDM!Mod(params.sefdm, SIGMA2);
    StopWatch swMod, swDem;

    while(totalbits < params.minTotalbits || errbits < params.maxErrbits )
    {
        if(totalbits > params.maxTotalbits)
            break;

        ++rndSeed;
        // auto chFR = makeChannel(params, rndSeed);

        // 送信ビット系列の生成
        Bit[] txbits;
        txbits = genBits(params.sefdm.nSpreadIn * params.sefdm.mod.symInputLength, rndSeed, txbits);

        // RDFT-s-SEFDM変調
        swMod.start();
        auto txsignal = modulateSEFDM(params.sefdm, txbits);
        swMod.stop();
        // auto rxsignal = txsignal.applyChannel(chFR);
        // rxsignal = rxsignal.addAWGN(rndSeed, SIGMA2);
        C[] rxsignal;
        params.channel.apply(txsignal, rxsignal);

        // OFDMによるチャネル推定
        // auto estChFR = estimateChannelByOFDM(params, chFR, SIGMA2);
        auto estChFR = estimateChannelByOFDM!C(params);

        swDem.start();
        Bit[] rxbits = demodulateSEFDM(params.sefdm, rxsignal, estChFR, params.channel.sigma2);
        // auto rxbits = sefdm.demodulate(txsignal);
        swDem.stop();
        errbits += txbits.zip(rxbits).map!"a[0] != a[1] ? 1 : 0".sum();
        totalbits += txbits.length;
    }


    SimResult result;
    result.totalbits = totalbits;
    result.errbits = errbits;
    result.ber = errbits * 1.0 / totalbits;
    result.modThroughput = totalbits * 1.0 / swMod.peek.total!"usecs" * 1e6;
    result.demodThroughput = totalbits * 1.0 / swDem.peek.total!"usecs" * 1e6;
    return result;
}


C[] modulateSEFDM(Mod)(in ref RDFTsSEFDMParams!Mod params, Bit[] txbits)
in(txbits.length == params.mod.symInputLength * params.nSpreadIn)
{
    // auto sefdm = new RDFTsSEFDM!(Mod, C)(params.mod, params.nFFT, 0, params.nTone, 1, params.nSpreadIn, params.nSpreadOut, null, params.rndSeed);
    auto sefdm = new RDFTsSEFDM!(Mod, C)(params.mod,
        params.nFFT, nCP: params.nFFT/4, params.nTone, nUpSampling: 1,
        params.nSpreadIn, params.nSpreadOut, shuffle: null, params.rndSeed, N0: 0, params.nIter);


    C[] dst;
    return sefdm.modulate(txbits, dst);
    // import std;
    // import dffdd.utils.fft;

    // C[] scmod;
    // scmod = params.mod.modulate(txbits, scmod);

    // auto fftw = makeFFTWObject!(TemplateOf!C)(params.nData);

    // // DFTスプレッド
    // fftw.inputs!F[] = scmod[];
    // fftw.fft!F();
    // scmod[] = fftw.outputs!F[];
    // foreach(ref e; scmod) {
    //     e /= sqrt(scmod.length * 1.0);          // FFTをユニタリにする
    //     e /= sqrt(params.nTone * 1.0 / params.nData);  // 1/sqrt(alpha)
    // }

    // // ランダムシャッフル
    // auto rnd = makeRNG(params.rndSeed, "RND_SHUF_SEFDM");
    // scmod = randomShuffle(scmod, rnd);

    // // 帯域圧縮
    // scmod.length = params.nFFT;
    // scmod[params.nTone .. $] = C(0);

    // fftw = makeFFTWObject!(TemplateOf!C)(params.nFFT);
    // fftw.inputs!F[] = scmod;
    // fftw.ifft!F();
    // scmod[] = fftw.outputs!F[];
    // foreach(ref e; scmod) {
    //     e *= sqrt(scmod.length * 1.0);          // IFFTをユニタリにする
    // }

    // return scmod;
}


Bit[] demodulateSEFDM(Mod)(in ref RDFTsSEFDMParams!Mod params, C[] rxsignal, C[] chFR, double N0)
{
    // auto sefdm = new RDFTsSEFDM!(Mod, C)(params.mod, params.nFFT, 0, params.nSpreadIn, params.nSpreadOut, 1, params.nData, null, params.rndSeed, N0, params.nIter);
    auto sefdm = new RDFTsSEFDM!(Mod, C)(params.mod,
        params.nFFT, nCP: params.nFFT/4, params.nTone, nUpSampling: 1,
        params.nSpreadIn, params.nSpreadOut, shuffle: null, params.rndSeed, N0, params.nIter);

    import mir.ndslice : sliced;
    import dffdd.math.vector : vectored;
    // sefdm.epDetector.channelSingularValues[] = chFR.sliced.vectored;
    sefdm.setFrequencyResponse(chFR);

    Bit[] dst;
    return sefdm.demodulate(rxsignal, dst);
}


Random makeRNG(T)(size_t seed, auto ref T id)
{
    import std;
    return Random(cast(uint) hashOf(forward!id, seed));
}


Bit[] genBits(size_t nbits, size_t seed, return ref Bit[] dst)
{
    import dffdd.utils.binary;
    import std;

    if(dst.length != nbits)
        dst.length = nbits;

    randomBits(makeRNG(seed, "GEN_BITS")).take(nbits).map!(a => Bit(a)).copy(dst);
    return dst;
}


C[] estimateChannelByOFDM(C, Mod)(ref SimParams!Mod params)
{
    immutable SIGMA2 = params.channel.sigma2;

    auto sefdm = new RDFTsSEFDM!(Mod, C)(params.sefdm.mod,
        params.sefdm.nFFT, nCP: params.sefdm.nFFT/4, params.sefdm.nTone, nUpSampling: 1,
        params.sefdm.nSpreadIn, params.sefdm.nSpreadOut, shuffle: null, params.sefdm.rndSeed, SIGMA2, params.sefdm.nIter);
    auto ofdm = sefdm.ofdm;

    // チャネル推定用のOFDM信号の生成
    C[] chEstSCtx = new C[params.sefdm.nTone];
    foreach(i; 0 .. params.sefdm.nTone)
        chEstSCtx[i] = i % 2 == 0 ? C(1) : C(-1);
    
    C[] chEstOFDMtx;
    ofdm.modulate(chEstSCtx, chEstOFDMtx);

    C[] estChFR = new C[params.sefdm.nTone];
    estChFR[] = C(0);
    foreach(i; 0 .. params.numSymChEst) {
        auto sig = chEstOFDMtx.dup;

        C[] rxsignal;
        if(params.chEstMethod == ChannelEstimationMethod.perfect) {
            // 雑音なしで完全に推定できると仮定される場合
            params.channel.apply(sig, rxsignal, Yes.ignoreAWGN);
        } else {
            // その他のときは，雑音を含めて受信された信号を使う
            params.channel.apply(sig, rxsignal);
        }

        C[] chEstSCrx;
        ofdm.demodulate(rxsignal, chEstSCrx);

        foreach(k; 0 .. params.sefdm.nTone) {
            if(params.chEstMethod == ChannelEstimationMethod.MMSE)
                estChFR[k] += chEstSCrx[k] / (chEstSCtx[k].sqAbs + SIGMA2) * chEstSCtx[k].conj;
            else if(params.chEstMethod == ChannelEstimationMethod.ZF || params.chEstMethod == ChannelEstimationMethod.perfect)
                estChFR[k] += chEstSCrx[k] / chEstSCtx[k];
            else {
                import std.exception;
                enforce(0, "Unknow channel estimation method");
            }
        }
    }


    foreach(k; 0 .. params.sefdm.nTone) {
        estChFR[k] /= params.numSymChEst;
    }

    return estChFR;
}


class WirelessChannel(C)
{
    import dffdd.filter.state : MultiFIRState;
    import std;
    import dffdd.blockdiagram.noise;


    this(in C[] taps, double sigma2, uint seed)
    {
        _rnd = BoxMuller!(Random, C)(makeRNG(seed, "GEN_AWGN"));
        _scale = sqrt(sigma2 / 2);

        _taps = taps.dup;
        _filt = MultiFIRState!C(1, taps.length);
        foreach(i; 0 .. taps.length)
            _filt.weight[i, 0] = taps[i];
    }


    void apply(in C[] inputs, ref C[] outputs, Flag!"ignoreAWGN" ignoreAWGN = No.ignoreAWGN)
    {
        if(outputs.length != inputs.length)
            outputs.length = inputs.length;

        foreach(i; 0 .. inputs.length) {
            _filt.update(inputs[i]);
            outputs[i] = _filt.output;

            if(!ignoreAWGN) {
                outputs[i] += _rnd.front * _scale;
                _rnd.popFront();
            }
        }
    }


    immutable(C)[] taps() const { return _taps; }
    double sigma2() const { return _scale^^2 * 2; }

  private:
    BoxMuller!(Random, C) _rnd;
    double _scale;
    immutable(C)[] _taps;
    MultiFIRState!C _filt;
}


WirelessChannel!C makeAWGNChannel(C)(double sigma2, uint seed)
{
    C[] taps = [C(1, 0)];
    return new WirelessChannel!C(taps, sigma2, seed);
}


WirelessChannel!C makeRayleighChannel(C)(size_t nTaps, double sigma2, uint seed)
{
    import std;
    import dffdd.blockdiagram.noise;

    auto rnd = makeRNG(seed, "CHANNEL");
    auto bm = BoxMuller!(Random, C)(rnd);

    C[] impResp = new C[nTaps];
    foreach(i; 0 .. nTaps) {
        impResp[i] = bm.front / sqrt(nTaps * 2.0);
        bm.popFront();
    }

    return new WirelessChannel!C(impResp, sigma2, seed);
}
