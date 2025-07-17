module dffdd.dsp.statistics;

import std.algorithm;
import std.array;
import std.numeric;
import dffdd.utils.fft;
import std.range;
import std.math;
import std.typecons;
import dffdd.blockdiagram.utils;
import dffdd.math.complex : RealPartType;


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


real calculateInBandSIC(R)(ref R r, real samplingFreq, size_t res, size_t nFFT, size_t nSubCarrier, size_t nOversampling, size_t avg = 32)
{
    alias C = typeof(ElementType!R.init[0]);
    C[] buf0 = new C[res * avg];
    C[] buf1 = new C[res * avg];

    foreach(i; 0 .. res * avg){
        auto front = r.front;
        buf0[i] = front[0];
        buf1[i] = front[1];
    }

    auto spec0 = calculatePowerSpectralDensity(buf0, samplingFreq, res, avg, null);
    auto spec1 = calculatePowerSpectralDensity(buf1, samplingFreq, res, avg, null);
    swapHalf(spec0);
    swapHalf(spec1);

    real p0 = 0,
         p1 = 0;
    foreach(i; 0 .. res / nOversampling * nSubCarrier / nFFT){
        p0 += spec0[i] + spec0[$-1-i];
        p1 += spec1[i] + spec1[$-1-i];
    }

    return p0 / p1;
}


struct FilterBankSpectrumAnalyzerImpl(C)
{
    import dffdd.filterbank;
    import dffdd.window;
    import dffdd.math.complex;

    alias R = RealPartType!C;

    this(size_t avg, size_t nCH, size_t nTap)
    {
        auto proto = designKaiserSincPulse!R(nCH, nTap, kaiserBetaFromStopbandAttdB(200), 1).map!(a => C(a, 0)).array();
        _afb = new DFTAnalysisFilterBank!C(nCH, proto);
        _nCH = nCH;
        _nTap = nTap;
        _delay = _afb.numOfDelay;
        _avg = avg;
        _buf1 = new C[_nCH];
        _buf2 = new C[_nCH];
        _psd = new R[_nCH];
        _psd[] = 0;
        _inst = makeInstrument!C(&this.consume);
    }


    void put(C e)
    {
        _inst.put(e);
    }


    void consume(FiberRange!C r)
    {
        foreach(i; 0 .. _delay / _nCH + _avg) {
            foreach(j; 0 .. _nCH) {
                assert(!r.empty);
                _buf1[j] = r.front;
                r.popFront();
            }

            _afb(_buf1, _buf2);

            if(i >= _delay / _nCH) {
                foreach(j, e; _buf2)
                    _psd[j] = (_psd[j] * _cnt + _buf2[j].sqAbs / _nCH) / (_cnt + 1);
                
                ++_cnt;
            }
        }
    }


    const(R)[] psd() const
    {
        return _psd;
    }


    size_t avgCount() const
    {
        return _cnt;
    }


  private:
    RefCounted!(FiberBuffer!C) _inst;
    DFTAnalysisFilterBank!C _afb;
    R[] _psd;
    C[] _buf1, _buf2;
    size_t _nCH, _nTap;
    size_t _avg;
    size_t _delay;
    size_t _cnt;
}


RefCounted!(FilterBankSpectrumAnalyzerImpl!C) makeFilterBankSpectrumAnalyzer(C)(size_t avg, size_t nCH, size_t nTap)
{
    return RefCounted!(FilterBankSpectrumAnalyzerImpl!C)(avg, nCH, nTap);
}
