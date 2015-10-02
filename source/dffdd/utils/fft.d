module dffdd.utils.fft;

import std.algorithm : swap;
import std.complex;
import std.math;
import std.numeric : Fft;
import std.stdio : File;

import carbon.math;

Complex!F[] fftWithSwap(F)(Fft fftObj, in Complex!F[] buf)
in{
    assert(buf.length.isPowOf2);
}
body{
    Complex!F[] spec = new Complex!F[buf.length];
    fftObj.fft(buf, spec);
    foreach(i; 0 .. spec.length/2)
        swap(spec[i], spec[$/2+i]);

    return spec;
}


C[] fftWithSwap(C)(Fft fftObj, in C[] buf)
if(is(typeof(buf[0].re) : float) && is(typeof(buf[0].im) : float))
in{
    assert(buf.length.isPowOf2);
}
body{
    alias F = typeof(buf[0].re);

    auto cpxBuf = new Complex!F[buf.length];
    foreach(i, ref e; cpxBuf)
        e = Complex!F(buf[i].re, buf[i].im);

    auto spec = fftObj.fftWithSwap(cpxBuf);
    C[] ret = new C[buf.length];

    foreach(i, ref e; ret)
        e = spec[i].re + spec[i].im * 1i;

    return ret;
}


align(1)
struct FrequencyDomain(C)
{
    C value;
    alias value this;
}


FrequencyDomain!C frequencyDomain(C)(C c)
{
    return typeof(return)(c);
}


inout(C)[] convAuto(C)(inout(FrequencyDomain!C)[] array)
{
    return (cast(inout(C)*)array.ptr)[0 .. array.length];
}

unittest
{
    FrequencyDomain!(Complex!float)[] carr;
    foreach(i; 0 .. 10) carr ~= complex!float(i, 0).frequencyDomain;

    auto carr2 = carr.convAuto!(Complex!float[]);
    foreach(i; 0 .. 10){
        assert(approxEqual(carr[i].re, carr2[i].re));
        assert(approxEqual(carr[i].im, carr2[i].im));
    }
}


cfloat[] rawReadComplex(File file, cfloat[] buf)
{
    return file.rawRead(buf);
}


Complex!float[] rawReadComplex(File file, cfloat[] buf, Complex!float[] output)
{
    auto res = file.rawReadComplex(buf);
    foreach(i, e; res)
        output[i] = complex!float(e.re, e.im);

    return output[0 .. res.length];
}


