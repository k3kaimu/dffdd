module dffdd.dsp.statistics;

import std.array;
import std.numeric;
import dffdd.utils.fft;


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
    cfloat[] buf = new cfloat[res];
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