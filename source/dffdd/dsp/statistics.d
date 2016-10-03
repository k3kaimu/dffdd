module dffdd.dsp.statistics;

import std.array;
import std.numeric;
import dffdd.utils.fft;
import std.range;
import std.math;


real[] calculatePowerSpectralDensity(R)(ref R r, real samplingFreq, size_t res, size_t avg = 32, real[] dstBuf = null)
{
    alias C = ElementType!R;

    C[] buf = new C[res];
    real[] fftRes;
    //fftRes[] = 0;
    if(dstBuf.length >= res)
        fftRes = dstBuf;
    else
        fftRes = new real[res];

    fftRes[] = 0;

    Fft fft = new Fft(res);

    foreach(n; 0 .. avg){
        foreach(i; 0 .. res){
            buf[i] = r.front;
            r.popFront();
        }

        auto spc = fft.fftWithSwap(buf);
        foreach(i; 0 .. res){
            fftRes[i] += spc[i].re^^2 + spc[i].im^^2;
        }
    }

    foreach(ref e; fftRes){
        e /= avg;
        e /= res;
        e /= samplingFreq;
    }


    return fftRes[0 .. res];
}


real calculateSIC(R)(ref R r, real samplingFreq, size_t res, size_t nFFT, size_t nSubCarrier, size_t nOversampling, size_t avg = 32)
{
    real sum1 = 0,
         sum2 = 0;

    foreach(i; 0 .. res * avg){
        auto front = r.front;
        sum1 += front[0].re^^2 + front[0].im^^2;
        sum2 += front[1].re^^2 + front[1].im^^2;

        r.popFront();
    }

    return sum1 / sum2;
}
