module impl;

import std.stdio;

import uhd.usrp;
import uhd.capi.types.metadata;
import carbon.channel;
import std.random;
import std.complex;
import std.concurrency;
import std.algorithm;
import std.range;
import core.thread;
import std.math;
import std.datetime;
import std.bitmanip;

import dffdd.blockdiagram.noise;
import dffdd.mod.ofdm;
import dffdd.mod.qpsk;
import dffdd.dsp.convolution;
import dffdd.utils.fft;

//import zmqd;
import preamble;
import spec;

immutable nCP = Constant.OFDM.nCP;
immutable nSC = Constant.OFDM.nSC;
immutable nFFT = Constant.OFDM.nFFT;
immutable nOS = Constant.nOverSampling;
immutable nSYM = (nCP + nFFT) * nOS;
immutable aScale = 5;

alias C = Complex!float;
alias ICA = immutable(C)[];


Preambles preambles;

static this()
{
    preambles = Preambles(false);
}



ICA makePilot()
{
    OFDM!C ofdm = new OFDM!C(Constant.OFDM.nFFT, Constant.OFDM.nCP, Constant.OFDM.nSC, Constant.nOverSampling);
    ICA inputs;
    Random rnd;
    rnd.seed(__LINE__);

    foreach(i; 0 .. Constant.OFDM.nSC)
        inputs ~= cast(C)std.complex.expi(uniform01(rnd) * 2 * PI);

    Complex!float[] dst;
    ofdm.modulate(inputs, dst);
    //foreach(ref e; dst)
    //    e *= aScale;

    return cast(ICA)dst;
    // return preambles.myEST;
}


ICA getPilot()
{
    static ICA pilot;

    if(pilot is null)
        pilot = makePilot();

    return pilot;
}


ICA getPreamble()
{
    return preambles.tgRTS;
}


C[] modulateImpl(Mod)(ref OFDM!C ofdm, ref Mod qpsk, const(ubyte)[] binary)
{
    immutable txbytes = binary.length;


    BitArray bitsArray = (){
        ubyte[] bins = binary.dup;
        while(bins.length % size_t.sizeof) bins ~= cast(ubyte)0;
        return BitArray(bins, txbytes * 8);
    }();

    ubyte[] bits;
    foreach(e; bitsArray)
        bits ~= e ? cast(ubyte)1 : cast(ubyte)0;


    C[] modQPSK;
    qpsk.modulate(bits, modQPSK);

    while(modQPSK.length % ofdm.symInputLength) modQPSK ~= C(0, 0);

    C[] dst;
    ofdm.modulate(modQPSK, dst);
    foreach(ref e; dst)
        e *= aScale;

    return dst;
}


C[] modulateSubcarriers(Mod)(ref OFDM!C ofdm, ref Mod qpsk, uint binSize, in Complex!float[] subcarriers)
{
    C[] txsignal;

    // プリアンブルを付加
    txsignal ~= getPreamble();

    // パイロット信号を追加
    auto pilot = getPilot.dup;
    foreach(ref e; pilot) e *= aScale;
    txsignal ~= pilot;

    // 最初にサイズ情報を付加
    uint len = binSize;
    txsignal ~= modulateImpl(ofdm, qpsk, (cast(ubyte*)&len)[0 .. uint.sizeof]);

    const(Complex!float)[] modQPSK = subcarriers;
    // ペイロードを付加
    while(modQPSK.length % ofdm.symInputLength) modQPSK ~= C(0, 0);

    C[] dst;
    ofdm.modulate(modQPSK, dst);
    foreach(ref e; dst)
        e *= aScale;

    return txsignal~dst;
}


C[] modulateTX(Mod)(ref OFDM!C ofdm, ref Mod qpsk, in ubyte[] binary, real payloadGain = 0.1)
{
    C[] txsignal;

    // プリアンブルを付加
    txsignal ~= getPreamble();

    // パイロット信号を追加
    txsignal ~= getPilot();

    // 最初にサイズ情報を付加
    uint len = cast(uint)binary.length;
    txsignal ~= modulateImpl(ofdm, qpsk, (cast(ubyte*)&len)[0 .. uint.sizeof]);

    // ペイロードを付加
    txsignal ~= (){
        auto s = modulateImpl(ofdm, qpsk, binary);
        immutable p = s.map!sqAbs.sum() / s.length;
        immutable g1 = sqrt(0.02 / p);
        immutable g2 = sqrt(payloadGain^^2 / p);
        foreach(ref e; s[0 .. $/2]) e *= g1;
        foreach(ref e; s[$/2 .. $]) e *= g2;
        return s;
    }();

    //txsignal ~= modulateImpl(ofdm, qpsk, binary);

    return txsignal;
}


size_t numOfHeaderSamples(Mod)(ref OFDM!C ofdm, ref Mod qpsk)
{
    immutable preamble = getPreamble().length;
    immutable pilot = getPilot().length;
    immutable sizeinfo = (){
        uint len = 0;
        return modulateImpl(ofdm, qpsk, (cast(ubyte*)&len)[0 .. uint.sizeof]).length;
    }();

    return preamble + pilot + sizeinfo;
}


size_t numOfPayloadSample(Mod)(ref OFDM!C ofdm, ref Mod qpsk, in ubyte[] binary)
{
    return modulateImpl(ofdm, qpsk, binary).length;
}


inout(C)[] findPreamble(PreambleDetector detector, inout(C)[] received)
{
    if(received.length < 4096 * 2){
        C[] buf = new C[4096 * 2];
        buf[] = C(0);
        buf[0 .. received.length] = received[];

        auto res = detector.detect!"TgRTS"(buf[0 .. $ / 4096 * 4096]);
        if(res[0] && res[1] < received.length)
            return received[res[1] .. $];
        else
            return null;
    }

    auto res = detector.detect!"TgRTS"(received[0 .. $ / 4096 * 4096]);

    if(!res[0])
        return null;
    else{
        writeln("Find Peak: ", res[1]);
        return received[res[1] .. $];
    }
}


unittest
{
    OFDM!C ofdm = new OFDM!C(Constant.OFDM.nFFT, Constant.OFDM.nCP, Constant.OFDM.nSC, Constant.nOverSampling);
    QPSK qpsk;
    PreambleDetector detector = new PreambleDetector(preambles);

    ubyte[] data;
    foreach(i; 0 .. 10)
        data ~= cast(ubyte)i;

    auto signal = modulateTX(ofdm, qpsk, data);
    foreach(ref e; signal)
         e *= C(0.5, -1.5);

    foreach(i; 1 .. signal.length)
        signal[i] = signal[i] + signal[i-1] * 0.01;

    size_t len = signal.length;

    C[] rnd;
    foreach(i; 0 .. 1025)
        rnd ~= C(uniform01, uniform01);

    signal = rnd ~ signal;
    assert(findPreamble(detector, signal).length == len);
}


BitArray demodulateImpl(Mod)(ref OFDM!C ofdm, ref Mod mod, in C[] freqResponse, in C[] symbols, ref C[] receivedSubcarriers)
{
    ofdm.demodulate(symbols, receivedSubcarriers);
    foreach(i; 0 .. receivedSubcarriers.length)
        receivedSubcarriers[i] = receivedSubcarriers[i] / freqResponse[i % $];

    //BitArray bits;
    ubyte[] bits = new ubyte[receivedSubcarriers.length / mod.symInputLength];
    mod.demodulate(receivedSubcarriers, bits);
    BitArray arr;
    foreach(e; bits)
        arr ~= cast(bool)e;

    return arr;
}


ubyte[] demodulateRX(Mod)(ref OFDM!C ofdm, ref Mod mod, const(C)[] received, ref C[] receivedSubcarriers)
{
    if(received.length < (nSYM * 2 + preambles.myEST.length))
        return null;

    received = received[preambles.myEST.length .. $];

    // パイロットを取り出す
    auto recvPilot = received[0 .. nSYM];
    C[] demodPilot;
    C[] demodPurePilot;
    ofdm.demodulate(recvPilot, demodPilot);
    ofdm.demodulate(getPilot, demodPurePilot);

    C[] freqResponse = new C[](nSC);
    foreach(i; 0 .. nSC)
        freqResponse[i] = demodPilot[i] / demodPurePilot[i];

    // バイト数が書かれている部分を取り出す
    C[] _dummy_;
    auto numBytes = (cast(uint[])cast(void[])demodulateImpl(ofdm, mod, freqResponse, received[nSYM .. nSYM * 2], _dummy_))[0];

    received = received[nSYM * 2 .. $];
    received = received[0 .. $ / nSYM * nSYM];

    auto bytes = (cast(ubyte[])cast(void[])demodulateImpl(ofdm, mod, freqResponse, received, receivedSubcarriers));
    if(bytes.length < numBytes)
        return null;
    else
        return bytes[0 .. numBytes].dup;
}


void demodulateSubcarriers(Mod)(ref OFDM!C ofdm, ref Mod mod, const(C)[] received, ref C[] receivedSubcarriers)
//if(is(Mod == BPSK))
{
    if(received.length < (nSYM * 2 + preambles.myEST.length))
        return;

    received = received[preambles.myEST.length .. $];

    // パイロットを取り出す
    auto recvPilot = received[0 .. nSYM];
    C[] demodPilot;
    C[] demodPurePilot;
    ofdm.demodulate(recvPilot, demodPilot);
    ofdm.demodulate(getPilot, demodPurePilot);

    C[] freqResponse = new C[](nSC);
    foreach(i; 0 .. nSC)
        freqResponse[i] = demodPilot[i] / demodPurePilot[i];

    // バイト数が書かれている部分を取り出す
    C[] _dummy_;
    auto numBytes = (cast(uint[])cast(void[])demodulateImpl(ofdm, mod, freqResponse, received[nSYM .. nSYM * 2], _dummy_))[0];

    received = received[nSYM * 2 .. $];
    received = received[0 .. $ / nSYM * nSYM];

    ofdm.demodulate(received, receivedSubcarriers);
    foreach(i; 0 .. receivedSubcarriers.length)
        receivedSubcarriers[i] = receivedSubcarriers[i] / freqResponse[i % $];

    if(receivedSubcarriers.length >= numBytes * 8)
        receivedSubcarriers = receivedSubcarriers[0 .. numBytes * 8];
}


unittest
{
    OFDM!C ofdm = new OFDM!C(Constant.OFDM.nFFT, Constant.OFDM.nCP, Constant.OFDM.nSC, Constant.nOverSampling);
    QPSK qpsk;

    ubyte[] data;
    foreach(i; 0 .. 111)
        data ~= cast(ubyte)i;

    auto signal = modulateTX(ofdm, qpsk, data);
    foreach(ref e; signal)
         e *= C(0.5, -1.5);

    foreach(i; 1 .. signal.length)
        signal[i] = signal[i] + signal[i-1];

    auto recvd = demodulateRX(ofdm, qpsk, signal);

    assert(recvd == data);
}
