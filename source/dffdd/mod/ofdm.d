module dffdd.mod.ofdm;

import carbon.math;
import std.stdio;

import dffdd.utils.fft;

final class OFDM(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;
    alias RealType = typeof(C.init.re);

    this(uint nFFT, uint nCp, uint nTone, uint nUpSampling = 1)
    {
        _nFFT = nFFT;
        _nCp = nCp;
        _nTone = nTone;
        _nUpSampling = nUpSampling;
        _inpBuffer = new C[nFFT * nUpSampling];
        _fftw = makeFFTWObject!C(nFFT * nUpSampling);
    }


    size_t symInputLength() const @property { return _nTone; }
    size_t symOutputLength() const @property { return (_nFFT + _nCp) * _nUpSampling; }


    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    in{
        assert(inputs.length % this.symInputLength == 0);
        assert(outputs.length % this.symOutputLength == 0);
        //assert(inputs.length / this.symInputLength == outputs.length / this.symOutputLength);
    }
    body{
        _inpBuffer[] = complexZero!C;

        if(outputs.length != inputs.length / this.symInputLength * this.symOutputLength)
            outputs.length = inputs.length / this.symInputLength * this.symOutputLength;

        auto mainLobeH = _inpBuffer[0 .. _nFFT/2];
        auto mainLobeL = _inpBuffer[$ - _nFFT/2 .. $];

        foreach(i; 0 .. inputs.length / _nTone){
            //_inpBuffer[] = complexZero!C;
            mainLobeL[] = complexZero!C;
            mainLobeH[] = complexZero!C;

            auto inpTones = inputs[i*_nTone .. (i+1)*_nTone];
            if(_nTone == _nFFT){
                mainLobeH[0 .. _nTone/2] = inpTones[0 .. _nTone/2];
                mainLobeL[$-_nTone/2 .. $] = inpTones[_nTone/2 .. $];

            }else if(_nTone == _nFFT -1){
                mainLobeH[1 .. $] = inpTones[0 .. _nFFT/2-1];
                mainLobeL[0 .. $] = inpTones[_nFFT/2-1 .. _nTone];
            }else if(_nTone % 2 == 0){
                mainLobeH[1 .. _nTone/2 + 1] = inpTones[0 .. _nTone/2];
                mainLobeL[$ - _nTone/2 .. $] = inpTones[_nTone/2 .. _nTone];
            }else{
                mainLobeH[1 .. (_nTone+1)/2+1] = inputs[0 .. (_nTone+1)/2];
                mainLobeL[$ - _nTone/2 .. $] = inpTones[(_nTone+1)/2 .. _nTone];
            }

            auto dst = outputs[i*symOutputLength .. (i+1)*symOutputLength];
            .ifft!RealType(_fftw, _inpBuffer, dst[_nUpSampling * _nCp .. _nUpSampling * (_nCp + _nFFT)]);
            dst[0 .. _nUpSampling * _nCp] = dst[$ - _nUpSampling * _nCp .. $];
            //dst[0] /= 4;
            //dst[1] /= 2;
            //dst[$-2] /= 2;
            //dst[$-1] /= 4;
            assert(dst.ptr == outputs.ptr + i*symOutputLength);
        }

        return outputs;
    }


    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    in{
        assert(outputs.length % this.symInputLength == 0);
        assert(inputs.length % this.symOutputLength == 0);
        //assert(outputs.length / this.symInputLength == inputs.length / this.symOutputLength);
    }
    body{
        //outputs[] = complexZero!C;

        auto mainLobeH = _inpBuffer[0 .. _nFFT/2];
        auto mainLobeL = _inpBuffer[$ - _nFFT/2 .. $];

        if(outputs.length != inputs.length / this.symOutputLength * this.symInputLength)
            outputs.length = inputs.length / this.symOutputLength * this.symInputLength;

        foreach(i; 0 .. inputs.length / this.symOutputLength){
            auto inputSymbol = inputs[i*this.symOutputLength .. (i+1)*this.symOutputLength];
            auto outputSymbol = outputs[i*this.symInputLength .. (i+1)*this.symInputLength];
            .fft!RealType(_fftw, inputSymbol[_nUpSampling * _nCp .. $], _inpBuffer);

            if(_nTone == _nFFT){
                outputSymbol[] = _inpBuffer[];
            }else if(_nTone == _nFFT -1){
                outputSymbol[0 .. _nFFT/2 -1] = mainLobeH[1 .. $];
                outputSymbol[_nFFT/2 -1.. _nTone] = mainLobeL[0 .. $];
            }else if(_nTone % 2 == 0){
                outputSymbol[0 .. _nTone/2] = mainLobeH[1 .. _nTone/2 +1];
                outputSymbol[_nTone/2 .. _nTone] = mainLobeL[$ - _nTone/2 .. $];
            }else{
                outputSymbol[0 .. (_nTone+1)/2] = mainLobeH[1 .. (_nTone+1)/2+1];
                outputSymbol[(_nTone+1)/2 .. _nTone] = mainLobeH[$ - _nTone/2 .. $];
            }
        }

        return outputs;
    }


  private:
    FFTWObject!C _fftw;
    C[] _inpBuffer;
    uint _nFFT;
    uint _nCp;
    uint _nTone;
    uint _nUpSampling;
}
