module dffdd.utils.fft;

import std.algorithm : swap;
import std.complex;
import std.math;
import std.numeric : Fft;

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
