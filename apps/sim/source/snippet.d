module snippet;


import std.stdio;
import std.math;
import std.complex;
import std.range;
import std.algorithm;

import dffdd.dsp.statistics;
import dranges.range;
import carbon.math;

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


auto connectToSwitch(R)(R r, const(bool)* sw, ElementType!R zero = complexZero!(ElementType!R))
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
        const(bool)* _sw;
        ElementType!R _zero;
    }


    return Result(r, sw, zero);
}


final class OFDMSymbolExchanger(Mod)
{
    alias InputElementType = Mod.InputElementType;
    alias OutputElementType = Mod.OutputElementType;


    this(Mod mod, const(bool)* sw, uint nFFT, uint nCp, uint nTone, uint nUpSampling = 1)
    {
        _mod = mod;
        _sw = sw;
        _nFFT = nFFT;
        _nCp = nCp;
        _nTone = nTone;
        _nUpSampling = nUpSampling;
    }


    size_t symInputLength() const @property { return _mod.symInputLength*2; }
    size_t symOutputLength() const @property { return _mod.symOutputLength*2; }


    ref OutputElementType[] modulate(in InputElementType[] input, return ref OutputElementType[] output)
    in{
        assert(input.length % this.symInputLength == 0);
    }
    body{
        _mod.modulate(input, output);

        if(*_sw){
            foreach(i; 0 .. input.length / this.symInputLength) {
                auto symAB = output[i * this.symOutputLength .. (i+1)*this.symOutputLength];

                auto startOfA = symAB[_nCp * _nUpSampling .. $];
                auto startOfB = symAB[$/2 + _nCp * _nUpSampling .. $];

                foreach(j; 0 .. _nFFT * _nUpSampling / 2)
                    swap(startOfA[j], startOfB[j]);
            }
        }

        return output;
    }


    ref InputElementType[] demodulate(in OutputElementType[] input, return ref InputElementType[] output)
    {
        if(output.length != input.length / this.symOutputLength * this.symInputLength)
            output.length = input.length / this.symOutputLength * this.symInputLength;

        auto input2 = input.dup;

        if(*_sw){
            foreach(i; 0 .. input.length / this.symOutputLength) {
                auto symAB = input2[i * this.symOutputLength .. (i+1)*this.symOutputLength];

                auto startOfA = symAB[_nCp * _nUpSampling .. $];
                auto startOfB = symAB[$/2 + _nCp * _nUpSampling .. $];

                foreach(j; 0 .. _nFFT * _nUpSampling / 2)
                    swap(startOfA[j], startOfB[j]);
            }
        }

        _mod.demodulate(input2, output);

        return output;
    }


  private:
    Mod _mod;
    const(bool)* _sw;
    uint _nFFT, _nCp, _nTone, _nUpSampling;
}


OFDMSymbolExchanger!Mod makeOFDMSymbolExchanger(Mod)(Mod mod, const(bool)* sw, uint nFFT, uint nCp, uint nTone, uint nUpSampling = 1)
{
    return new OFDMSymbolExchanger!Mod(mod, sw, nFFT, nCp, nTone, nUpSampling);
}