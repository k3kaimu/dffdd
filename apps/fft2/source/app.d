import std.stdio;
import std.stdio;
import std.getopt;
import std.range;
import std.algorithm;
import std.array;
import std.format;
import std.string;
import std.numeric;
import std.math;
import std.complex;
import std.datetime;

import carbon.stream;
import carbon.math;

import dffdd.filter.lms;
import dffdd.filter.polynomial;
import dffdd.filter.diagonal;
import dffdd.filter.mempoly;

import dffdd.utils.fft;
enum real sampFreq = 400e3;


void main(string[] args)
{
    auto file = File(args[1]);

    file.seek(919 * 8);
    enum size_t BufSize = (cast(size_t)sampFreq).nextPowOf2;

    cfloat[] buf = new cfloat[BufSize];
    file.rawRead(buf);
    foreach(ref e; buf)
        e *= e;

    Fft fft = new Fft(BufSize);
    auto fftRes = fft.fftWithSwap(buf);

    float maxP = 0; size_t maxIndex;
    foreach(i, e; fftRes)
    {
        float p = e.re^^2 + e.im^^2;
        if(maxP < p){maxP = p; maxIndex = i;}
    }

    writefln("%s: %s[Hz]", maxP, maxIndex*sampFreq/BufSize-(sampFreq/2));
}
