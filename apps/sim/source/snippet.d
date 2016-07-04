module snippet;


import std.stdio;
import std.math;
import std.complex;
import std.range;

import dffdd.dsp.statistics;
import dranges.range;

void writePSD(R)(auto ref R r, File file, real samplingFreq, size_t resolution)
{
    auto psd = r.calculatePowerSpectralDensity(samplingFreq, resolution);
    foreach(i, e; psd){
        file.writefln("%f,%f,%f,", (i*1.0/resolution*samplingFreq-(samplingFreq/2))/1e6, e, e == 0 ? -400 : 30+10*log10(e));
    }
}


real measurePower(R)(auto ref R r, size_t n)
{
    real s = 0;
    foreach(i; 0 .. n){
        if(r.empty){
            n = i;
            break;
        }

        auto e = r.front;
        s += e.re^^2 + e.im^^2;
        r.popFront();
    }

    return s / n;
}


real measureBER(R1, R2)(ref R1 r1, ref R2 r2, ulong totalBits)
{
    r1.popFrontN(10*1000);
    r2.popFrontN(10*1000);

    ulong cnt;
    foreach(i; 0 .. totalBits){
        auto a = r1.front;
        auto b = r2.front;

        if(a != b) ++cnt;

        r1.popFront();
        r2.popFront();
    }

    return (cast(real)cnt)/totalBits;
}


auto connectToSwitch(R)(R r, bool* sw, ElementType!R zero = 0+0i)
{
    static struct Result 
    {
        auto front() @property
        {
            if(*_sw)
                return _r.front;
            else
                return _zero;
        }


        void popFront()
        {
            if(*_sw) _r.popFront();
        }


        bool empty()
        {
            return *_sw && _r.empty;
        }

      private:
        R _r;
        bool* _sw;
        ElementType!R _zero;
    }


    return Result(r, sw, zero);
}
