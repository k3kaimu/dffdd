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
//enum real numOfSumple = 1024;
import dffdd.utils.unit;
import dffdd.blockdiagram.amplifier;
import dffdd.blockdiagram.iqmixer;
import dffdd.blockdiagram.noise;
import dffdd.blockdiagram.utils;
import dffdd.blockdiagram.vga;
import dffdd.blockdiagram.adder;
import dffdd.blockdiagram.decimator;
import dffdd.blockdiagram.quantizer;
import dffdd.gps.code;


enum real samplingFreq = 40e6;
enum size_t upSamplingRatio = 1;
enum real upSampledFreq = upSamplingRatio * samplingFreq;
enum real ncoFreq = 8e6;
enum size_t Res = 1024;


void main()
{
    auto code1 = L1CACode(4),
         code2 = L1CACode(5);

    writeln(code1.zip(code2).map!"a[0]*a[1]".take(1023).sum());

/+
    auto y = ThermalNoise(SamplingFrequency(samplingFreq), 300);
    //auto y = boxMullerNoise().map!"a.re+0i";
    //auto y = repeat(1.0+0*1i);
    //auto y = preciseComplexNCO(1e6/*+ncoFreq/4*/, 1.0/upSampledFreq, 0)
    //.add(ThermalNoise(SamplingFrequency(samplingFreq), 300))
    //.add(preciseComplexNCO(10e6+ncoFreq/8, 1.0/upSampledFreq, 0));
    //.connectTo!VGA((-50).dB)
    //.connectTo!IQImbalance(0.dB, 25.dB)
    //.connectTo!PhaseNoise((-30).dB, 0.1)
    //.connectTo!VGA(10.dB)
    //.connectTo!PowerAmplifier(30.dB, 10.dBm)
    //.connectToPowerOf2Decimator!upSamplingRatio(1024)
    //.connectTo!Quantizer(4);
    //.connectTo!(ADC!16)();

    foreach(i; 0 .. Res)
        y.popFront();

    cfloat[] buf = new cfloat[Res];
    real[] fftRes = new real[Res];
    fftRes[] = 0;

    foreach(n; 0 .. 32){
        foreach(i; 0 .. Res){
            buf[i] = y.front;
            y.popFront();
        }

        Fft fft = new Fft(Res);

        auto res = fft.fftWithSwap(buf);
        foreach(i; 0 .. Res){
            fftRes[i] += res[i].re^^2 + res[i].im^^2;
        }
    }

    foreach(ref e; fftRes){
        e /= 32;
        e /= samplingFreq;
        e /= Res;
    }

    size_t i = 0;
    foreach(e; fftRes)
    {
        writefln("%f,%f,%f", (i*1.0/Res*samplingFreq-(samplingFreq/2))/1e6, e, e == 0 ? -400 : 30+10*log10(e));
        ++i;
    }

    //auto file = File(args[1]);

    //file.seek(100 * 8);
    //enum size_t BufSize = (cast(size_t)numOfSumple).nextPowOf2;

    //cfloat[] buf = new cfloat[BufSize];
    //file.rawRead(buf);

    //Fft fft = new Fft(BufSize);
    //auto fftRes = fft.fftWithSwap(buf);

    //foreach(e; fftRes.map!"a.re^^2 + a.re^^2")
    //{
    //    writefln("%f,%f", e, 10*log10(e));
    //}
    +/
}
