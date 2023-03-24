module app;

import std.algorithm;
import std.complex;
import std.random;
import std.range;
import std.stdio;
import std.traits;

import dffdd.mod;


alias F = double;
alias C = Complex!F;


void main()
{
    // writeln("Hello");
    foreach(vEbN0dB; [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) {
        SimParams!(QPSK!C) params;
        params.sefdm.nFFT = 64*4*4;
        params.sefdm.nTone = 40*4*4;
        params.sefdm.rndSeed = 0;
        params.vEbN0dB = vEbN0dB;

        auto simResult = runSimImpl(params);
        writefln!"%s: %s, Mod: %s [Mbps], Demod: %s [Mbps]"(vEbN0dB, simResult.ber, simResult.modThroughput / 1e6, simResult.demodThroughput / 1e6);
    }
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
    size_t nFFT;
    size_t nTone;
    size_t rndSeed = 0;

    alias nData = nFFT;
}


struct SimParams(Mod)
{
    RDFTsSEFDMParams!Mod sefdm;
    double vEbN0dB;
    size_t maxErrbits = 100;
    size_t minTotalbits = 100_000_000;
}


SimResult runSimImpl(Mod)(SimParams!Mod params)
{
    import std.datetime.stopwatch;

    immutable double SIGMA2 = 1.0 / (10^^(params.vEbN0dB/10.0) * params.sefdm.mod.symInputLength);

    size_t totalbits = 0;
    size_t errbits = 0;

    size_t rndSeed = 0;

    // auto sefdm = new RDFTsSEFDM!Mod(params.sefdm, SIGMA2);
    StopWatch swMod, swDem;

    while(totalbits < params.minTotalbits && errbits < params.maxErrbits)
    {
        ++rndSeed;

        // 送信ビット系列の生成
        swMod.start();
        Bit[] txbits;
        txbits = genBits(params.sefdm.nData * params.sefdm.mod.symInputLength, rndSeed, txbits);
        swMod.stop();

        // RDFT-s-SEFDM変調
        auto txsignal = modulateSEFDM(params.sefdm, txbits);
        // auto txsignal = sefdm.modulate(txbits);
        txsignal = txsignal.addAWGN(rndSeed, SIGMA2);

        swDem.start();
        Bit[] rxbits = demodulateSEFDM(params.sefdm, txsignal, SIGMA2);
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
in(txbits.length == params.mod.symInputLength * params.nData)
{
    import std;
    import dffdd.utils.fft;

    C[] scmod;
    scmod = params.mod.modulate(txbits, scmod);

    // DFTスプレッド
    auto fftw = makeFFTWObject!(TemplateOf!C)(params.nData);
    fftw.inputs!F[] = scmod[];
    fftw.fft!F();
    scmod[] = fftw.outputs!F[];
    foreach(ref e; scmod) {
        e /= sqrt(scmod.length * 1.0);          // FFTをユニタリにする
        e /= sqrt(params.nTone * 1.0 / params.nData);  // 1/sqrt(alpha)
    }

    // ランダムシャッフル
    auto rnd = makeRNG(params.rndSeed, "RND_SHUF_SEFDM");
    scmod = randomShuffle(scmod, rnd);

    // 帯域圧縮
    scmod.length = params.nFFT;
    scmod[params.nTone .. $] = C(0);

    fftw = makeFFTWObject!(TemplateOf!C)(params.nFFT);
    fftw.inputs!F[] = scmod;
    fftw.ifft!F();
    scmod[] = fftw.outputs!F[];
    foreach(ref e; scmod) {
        e *= sqrt(scmod.length * 1.0);          // IFFTをユニタリにする
    }

    return scmod;
}


Bit[] demodulateSEFDM(Mod)(in ref RDFTsSEFDMParams!Mod params, C[] rxsignal, double N0)
{
    import std;
    import mir.ndslice : sliced;
    import dffdd.math;
    import dffdd.detector.ep;
    import dffdd.utils.fft;

    // DFTスプレッド
    auto W1 = dftMatrix!C(params.nData);

    // ランダムシャッフル
    auto rnd = makeRNG(params.rndSeed, "RND_SHUF_SEFDM");
    size_t[] plist = iota(params.nData).map!"cast(size_t)a".array();
    plist = randomShuffle(plist, rnd);
    auto P = PermutationMatrix(plist);

    // 圧縮&電力割当
    C[] dlist = new C[params.nTone];
    dlist[] = C(1/sqrt(params.nTone * 1.0 / params.nData));

    // auto W2 = idftMatrix!C(params.nFFT);
    // EP復調器
    auto detector = makeSVDEPDetector!C(params.mod, identity!C(params.nTone), dlist.sliced.vectored, P * W1, N0, 20);

    // 受信信号を周波数領域へ変換
    {
        auto fftw = makeFFTWObject!(TemplateOf!C)(params.nFFT);
        fftw.inputs!F[] = rxsignal[];
        fftw.fft!F();
        rxsignal[] = fftw.outputs!F[];
        foreach(ref e; rxsignal)
            e /= sqrt(params.nFFT * 1.0);
    }

    C[] detected;
    detected = detector.detect(rxsignal[0 .. params.nTone], detected);

    Bit[] demodBits;
    demodBits = params.mod.demodulate(detected, demodBits);
    return demodBits;
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


R addAWGN(R)(R received, size_t seed, double SIGMA)
{
    import std;
    import dffdd.blockdiagram.noise;
    auto rnd = BoxMuller!(Random, C)(makeRNG(seed, "GEN_AWGN"));

    foreach(ref e; received) {
        e += rnd.front * sqrt(SIGMA/2);
        rnd.popFront();
    }

    return received;
}