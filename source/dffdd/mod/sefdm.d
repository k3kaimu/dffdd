module dffdd.mod.sefdm;

__EOF__

import std.random;
import std.traits;

import dffdd.mod.ofdm;
import dffdd.math.complex;
import dffdd.utils.fft;



final class RDFTsSEFDM(Mod, C = Mod.OutputELementType)
{
    alias InputElementType = Bit;
    alias OutputElementType = C;
    alias RealType = typeof(C.init.re);


    this(Mod mod, uint nFFT, uint nCP, uint nSC, uint nTone, uint nData, uint randSeed, uint nUpSampling = 1)
    {
        _mod = mod;
        _ofdm = new OFDM!C(nFFT, nCP, nSC, nTone, nUpSampling);
        _nTone = nTone;
        _nData = nData;
        _cbuffer = new C[_nData];
        _sbuffer = new C[_nTone];

        _fftw = makeFFTWObject!C(_nData);

        // シャッフルするインデックス配列
        _rndidx = new size_t[nData];
        foreach(i; 0 .. nData)
            _rndIdx[i] = i;
        
        Random rnd;
        rnd.seed(randSeed);
        _rndIdx = randomShuffle(_rndIdx, rnd);
    }


    size_t symInputLength() const @property { return _nData * mod.symInputLength; }
    size_t symOutputLength() const @property { return _ofdm.symOutputLength; }


    auto makeDetector(Mod)(Mod mod, in C[] chFreqResp, RealType sigma2)
    {
        auto matU = identity!C()
    }


    ref C[] modulate(in Bit[] inputs, return ref C[] outputs)
    in(inputs.length % this.symInputLength == 0)
    do {
        if(outputs.length != inputs.length / this.symInputLength * this.symOutputLength)
            outputs.length = inputs.length / this.symInputLength * this.symOutputLength;

        immutable nIN = this.symInputLength,
                  nOUT = this.symOutputLength;

        foreach(i; 0 .. inputs.length / _nData) {
            auto inbuf = inputs[i * nIN .. (i+1) * nIN];
            auto outbuf = outputs[i * nOUT .. (i+1) * nOUT];

            _cbuffer = _mod.modulate(inbuf, _cbuffer);
            _sbuffer = _randomSpread(_cbuffer, _sbuffer);
            _ofdm.modulate(_sbuffer, outbuf);
        }

        return outputs;
    }


    ref C[] demodulateAsOFDM(in C[] received, return ref C[] subcarriers)
    in(received.length % this.symOutputLength == 0)
    {
        return _ofdm.demodulate(received, subcarriers);      
    }


  private:
    alias F = typeof(C.init.re);

    OFDM!C _ofdm;
    FFTWObject!(TemplateOf!C) _fftw;
    C[] _sbuffer;
    size_t _nTone;
    size_t _nData;
    size_t[] _rndidx;


    void _randomSpread(in C[] inbuf, C[] outbuf)
    in(inbuf.length == _nData)
    in(outbuf.length == _nTone);
    {
        alias F = typeof(C.init.re);
        immutable F scale = 1/fast_sqrt!F(inbuf.length);

        _fftw.inputs!F[] = inbuf[];
        _fftw.fft!F();
        foreach(i; 0 .. _nTone)
            outbuf[i] = _fftw.outputs!F[_rndIdx[i]] * scale;
    }
}
