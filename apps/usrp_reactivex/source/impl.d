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

import dffdd.mod.ofdm;
import dffdd.mod.qpsk;
import dffdd.dsp.convolution;
import dffdd.utils.fft;

import zmqd;
import preamble;
import spec;

immutable nCP = 16;
immutable nSC = 52;
immutable nFFT = 64;
immutable nOS = 4;
immutable nSYM = (nCP + nFFT) * nOS;

alias C = Complex!float;
alias ICA = immutable(C)[];

Preambles preambles;

static this()
{
    preambles = Preambles(false);
}



ICA makePilot()
{
    OFDM!C ofdm = new OFDM!C(64, 16, 52, 4);
    ICA inputs;
    Random rnd;
    rnd.seed(__LINE__);

    foreach(i; 0 .. 52)
        inputs ~= cast(C)std.complex.expi(uniform01(rnd) * 2 * PI);

    Complex!float[] dst;
    ofdm.modulate(inputs, dst);

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


C[] modulateImpl(ref OFDM!C ofdm, ref QPSK qpsk, const(ubyte)[] binary)
{
    immutable txbytes = binary.length;


    BitArray bits = (){
        ubyte[] bins = binary.dup;
        while(bins.length % size_t.sizeof) bins ~= cast(ubyte)0;
        return BitArray(bins, txbytes * 8);
    }();


    C[] modQPSK;
    qpsk.modulate(bits, modQPSK);

    while(modQPSK.length % ofdm.symInputLength) modQPSK ~= C(0, 0);

    C[] dst;
    ofdm.modulate(modQPSK, dst);

    return dst;
}


C[] modulateTX(ref OFDM!C ofdm, ref QPSK qpsk, in ubyte[] binary)
{
    C[] txsignal;

    // プリアンブルを付加
    txsignal ~= getPreamble();

    // パイロット信号を追加
    txsignal ~= getPilot();

    // 最初にサイズ情報を付加
    size_t len = binary.length;
    txsignal ~= modulateImpl(ofdm, qpsk, (cast(ubyte*)&len)[0 .. size_t.sizeof]);

    // ペイロードを付加
    txsignal ~= modulateImpl(ofdm, qpsk, binary);

    return txsignal;
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
    else
        return received[res[1] .. $];
}


unittest
{
    OFDM!C ofdm = new OFDM!C(64, 16, 52, 4);
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


BitArray demodulateImpl(ref OFDM!C ofdm, ref QPSK qpsk, in C[] freqResponse, in C[] symbols)
{
    C[] subcarriers;
    ofdm.demodulate(symbols, subcarriers);
    foreach(i; 0 .. subcarriers.length)
        subcarriers[i] = subcarriers[i] / freqResponse[i % $];

    BitArray bits;
    qpsk.demodulate(subcarriers, bits);
    return bits;
}


ubyte[] demodulateRX(ref OFDM!C ofdm, ref QPSK qpsk, const(C)[] received)
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
    auto numBytes = (cast(size_t[])cast(void[])demodulateImpl(ofdm, qpsk, freqResponse, received[nSYM .. nSYM * 2]))[0];

    received = received[nSYM * 2 .. $];
    received = received[0 .. $ / nSYM * nSYM];

    return (cast(ubyte[])cast(void[])demodulateImpl(ofdm, qpsk, freqResponse, received))[0 .. numBytes].dup;
}


unittest
{
    OFDM!C ofdm = new OFDM!C(64, 16, 52, 4);
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
