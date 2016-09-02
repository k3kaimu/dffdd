module dffdd.dsp.statistics;

import std.array;
import std.numeric;
import dffdd.utils.fft;
import std.range;
import std.math;


//enum 

/+
auto average(alias add = "a+b", R)(ref R r)
if(!isInfinite!R)
{
    if(!r.empty){
        auto fd = r.fold!(add, "a+1")();
        return r[0] / r[1];
    }else
        return typeof(return).init;
}


auto average(alias add = "a+b", R)(ref R r, size_t len = 1024)
if(isInfinite!R)
{
    return r.take(1024).fold!add() / len;
}+/





//auto averager(HowToAverage)


//float blockAvera

/+
float averagePower(R)(ref R r, size_t len = 1024)
{
    float p = 0;
    foreach(i; 0 .. len){
        auto e = r.front;
        p += e.re^^2 + e.im^^2;
        r.popFront();
    }

    return (p / len);
}


auto blockAveragePower(R)(R r, size_t len = 1024)
{
    return r.chunks(len).map!(a => a.averagePower(len));
}


auto movingAveragePower(R)(R r, size_t len = 1024)
+/

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
    alias C = typeof(ElementType!R.init.tupleof[0]);

    C[] buf1 = new C[res];
    C[] buf2 = new C[res];

    real[] fftRes1 = new real[res],
           fftRes2 = new real[res];

    fftRes1[] = 0;
    fftRes2[] = 0;

    Fft fft = new Fft(res);

    foreach(n; 0 .. avg){
        foreach(i; 0 .. res){
            auto e = r.front;
            buf1[i] = e[0];
            buf2[i] = e[1];
            r.popFront();
        }

        auto spc1 = fft.fftWithSwap(buf1);
        auto spc2 = fft.fftWithSwap(buf2);

        foreach(i; 0 .. res){
            fftRes1[i] += spc1[i].re^^2 + spc1[i].im^^2;
            fftRes2[i] += spc2[i].re^^2 + spc2[i].im^^2;
        }
    }

    immutable boundF1 = samplingFreq / nOversampling / nFFT * nSubCarrier / 2 * 0.8;
    immutable boundF2 = samplingFreq / nOversampling / nFFT * 2;

    real sum1 = 0,
         sum2 = 0;
    foreach(i; 0 .. res){
        immutable f = abs(i * samplingFreq / res - samplingFreq/2);
        if(f < boundF1 && f > boundF2){
            sum1 += fftRes1[i];
            sum2 += fftRes2[i];
        }
    }

    return sum1 / sum2;
}
