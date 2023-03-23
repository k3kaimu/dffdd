module app;

import std.algorithm;
import std.complex;
import std.random;
import std.range;
import std.stdio;

import dffdd.mod;


alias F = double;
alias C = Complex!F;


void main()
{
    // writeln("Hello");
    foreach(vEbN0dB; [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) {
        SimParams!(QPSK!C) params;
        params.sefdm.nFFT = 64;
        params.sefdm.nTone = 64;
        params.sefdm.rndSeed = 0;
        params.vEbN0dB = vEbN0dB;

        auto simResult = runSimImpl(params);
        writefln!"%s: %s"(vEbN0dB, simResult.ber);
    }
}


struct SimResult
{
    double ber;
    long errbits;
    long totalbits;
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
    size_t minTotalbits = 10_000_000;
}


SimResult runSimImpl(Mod)(SimParams!Mod params)
{
    immutable double SIGMA2 = 1.0 / (10^^(params.vEbN0dB/10.0) * params.sefdm.mod.symInputLength);

    size_t totalbits = 0;
    size_t errbits = 0;

    size_t rndSeed = 0;

    while(totalbits < params.minTotalbits && errbits < params.maxErrbits)
    {
        ++rndSeed;

        // 送信ビット系列の生成
        Bit[] txbits;
        txbits = genBits(params.sefdm.nData * params.sefdm.mod.symInputLength, rndSeed, txbits);

        // RDFT-s-SEFDM変調
        auto txsignal = modulateSEFDM(params.sefdm, txbits);
        txsignal = txsignal.addAWGN(rndSeed, SIGMA2);

        Bit[] rxbits = demodulateSEFDM(params.sefdm, txsignal);
        errbits += txbits.zip(rxbits).map!"a[0] != a[1] ? 1 : 0".sum();
        totalbits += txbits.length;
    }


    SimResult result;
    result.totalbits = totalbits;
    result.errbits = errbits;
    result.ber = errbits * 1.0 / totalbits;
    return result;
}


C[] modulateSEFDM(Mod)(in ref RDFTsSEFDMParams!Mod params, Bit[] txbits)
in(txbits.length == params.mod.symInputLength * params.nData)
{
    C[] dst;
    return params.mod.modulate(txbits, dst);
}


Bit[] demodulateSEFDM(Mod)(in ref RDFTsSEFDMParams!Mod params, C[] rxsignal)
{
    Bit[] dst;
    return params.mod.demodulate(rxsignal, dst);
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