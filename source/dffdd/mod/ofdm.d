module dffdd.mod.ofdm;

import carbon.math;

import dffdd.utils.fft;

import std.math;
import dffdd.math : isNarrowComplex;

final class OFDM(C)
if(isNarrowComplex!C)
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
        _scaleTX = sqrt((nFFT * nUpSampling * 1.0)^^2 / nTone);
    }


    size_t symInputLength() const @property { return _nTone; }
    size_t symOutputLength() const @property { return (_nFFT + _nCp) * _nUpSampling; }


    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    in{
        assert(inputs.length % this.symInputLength == 0);
    }
    do{
        _inpBuffer[] = C(0);

        if(outputs.length != inputs.length / this.symInputLength * this.symOutputLength)
            outputs.length = inputs.length / this.symInputLength * this.symOutputLength;

        auto mainLobeH = _inpBuffer[0 .. _nFFT/2];
        auto mainLobeL = _inpBuffer[$ - _nFFT/2 .. $];

        foreach(i; 0 .. inputs.length / _nTone){
            mainLobeL[] = C(0);
            mainLobeH[] = C(0);

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
                mainLobeH[1 .. (_nTone+1)/2+1] = inpTones[0 .. (_nTone+1)/2];
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

        foreach(ref e; outputs)
            e *= _scaleTX;

        return outputs;
    }


    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    in{
        assert(inputs.length % this.symOutputLength == 0);
    }
    do{
        //outputs[] = C(0);

        auto mainLobeH = _inpBuffer[0 .. _nFFT/2];
        auto mainLobeL = _inpBuffer[$ - _nFFT/2 .. $];

        if(outputs.length != inputs.length / this.symOutputLength * this.symInputLength)
            outputs.length = inputs.length / this.symOutputLength * this.symInputLength;

        foreach(i; 0 .. inputs.length / this.symOutputLength){
            auto inputSymbol = inputs[i*this.symOutputLength .. (i+1)*this.symOutputLength];
            auto outputSymbol = outputs[i*this.symInputLength .. (i+1)*this.symInputLength];
            .fft!RealType(_fftw, inputSymbol[_nUpSampling * _nCp .. $], _inpBuffer);

            if(_nTone == _nFFT){
                outputSymbol[0 .. _nFFT/2] = mainLobeH[0 .. $];
                outputSymbol[_nFFT/2 .. _nFFT] = mainLobeL[0 .. $];
            }else if(_nTone == _nFFT -1){
                outputSymbol[0 .. _nFFT/2 -1] = mainLobeH[1 .. $];
                outputSymbol[_nFFT/2 -1.. _nTone] = mainLobeL[0 .. $];
            }else if(_nTone % 2 == 0){
                outputSymbol[0 .. _nTone/2] = mainLobeH[1 .. _nTone/2 +1];
                outputSymbol[_nTone/2 .. _nTone] = mainLobeL[$ - _nTone/2 .. $];
            }else{
                outputSymbol[0 .. (_nTone+1)/2] = mainLobeH[1 .. (_nTone+1)/2+1];
                outputSymbol[(_nTone+1)/2 .. _nTone] = mainLobeL[$ - _nTone/2 .. $];
            }
        }

        foreach(ref e; outputs)
            e /= _scaleTX;

        return outputs;
    }


  private:
    FFTWObject!C _fftw;
    C[] _inpBuffer;
    uint _nFFT;
    uint _nCp;
    uint _nTone;
    uint _nUpSampling;
    real _scaleTX;
}

//
unittest
{
    import std.complex, std.math;

    // FFTサイズ: 8
    // サイクリックプレフィックス: 3
    // 使用サブキャリア数: 4
    // アップサンプリング率: 2
    auto ofdmMod = new OFDM!(Complex!float)(8, 3, 4, 2);

    Complex!float[] inps = [Complex!float(1, 1), Complex!float(-1, -1), Complex!float(1, 0), Complex!float(0, 1)];
    Complex!float[] res;

    // 変調
    ofdmMod.modulate(inps, res);

    assert(res.length == (8 + 3) * 2);

    // CPのチェック
    foreach(i; 0 .. 3 * 2)
        assert(res[i] == res[$ - 3*2 + i]);

    Complex!float[] demodsym;
    ofdmMod.demodulate(res, demodsym);

    assert(demodsym.length == inps.length);

    foreach(i; 0 .. inps.length){
        assert(isClose(inps[i].re, demodsym[i].re));
        assert(isClose(inps[i].im, demodsym[i].im));
    }
}



/**
指定された位置にパイロットサブキャリアを挿入するOFDM変復調
ただし，パイロットサブキャリアを挿入する位置は固定で，サブキャリアの最も外側に1つずつ計2本配置します．
*/
final class OFDMWithStaticPilot(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;
    alias RealType = typeof(C.init.re);


    this(uint nFFT, uint nCp, uint nDataTone, uint nUpSampling)
    in {
        assert(nDataTone % 2 == 0);
    }
    do {
        _nDataTone = nDataTone;
        _ofdmMod = new OFDM!C(nFFT, nCp, nDataTone + 2, nUpSampling);

        _inpBuffer.length = nDataTone + 2;
    }


    size_t symInputLength() const @property { return _nDataTone; }
    size_t symOutputLength() const @property { return _ofdmMod.symOutputLength; }


    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    in {
        assert(inputs.length % this.symInputLength == 0);
    }
    do {
        if(outputs.length != inputs.length / this.symInputLength * this.symOutputLength)
            outputs.length = inputs.length / this.symInputLength * this.symOutputLength;

        _inpBuffer[$/2 - 1] = 1;
        _inpBuffer[$/2] = -1;
        auto mainLobeH = _inpBuffer[0 .. $/2 - 1];
        auto mainLobeL = _inpBuffer[$/2 + 1 .. $];

        foreach(i; 0 .. inputs.length / _nDataTone){
            auto inpTones = inputs[i*_nDataTone .. (i+1)*_nDataTone];
            auto outSignal = outputs[i*this.symOutputLength .. (i+1)*this.symOutputLength];

            mainLobeH[] = inpTones[0 .. $/2];
            mainLobeL[] = inpTones[$/2 .. $];

            _ofdmMod.modulate(_inpBuffer, outSignal);
        }

        return outputs;
    }


    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    in {
        assert(inputs.length % this.symOutputLength == 0);
    }
    do {
        auto mainLobeH = _inpBuffer[0 .. $/2 - 1];
        auto mainLobeL = _inpBuffer[$/2 + 1 .. $];

        if(outputs.length != inputs.length / this.symOutputLength * this.symInputLength)
            outputs.length = inputs.length / this.symOutputLength * this.symInputLength;

        foreach(i; 0 .. inputs.length / this.symOutputLength){
            auto inputSymbol = inputs[i*this.symOutputLength .. (i+1)*this.symOutputLength];
            auto outputSymbol = outputs[i*this.symInputLength .. (i+1)*this.symInputLength];

            _ofdmMod.demodulate(inputSymbol, _inpBuffer);
            C coef = (_inpBuffer[$/2-1] - _inpBuffer[$/2])/2;
            foreach(ref e; _inpBuffer) e /= coef;

            outputSymbol[0 .. $/2] = mainLobeH[];
            outputSymbol[$/2 .. $] = mainLobeL[];
        }

        return outputs;
    }


  private:
    size_t _nDataTone;
    immutable(int)[] _pilotIndex;
    OFDM!C _ofdmMod;
    C[] _inpBuffer;
}


//
unittest
{
    import std.complex, std.math;

    // FFTサイズ: 8
    // サイクリックプレフィックス: 3
    // 使用データサブキャリア数: 4
    // アップサンプリング率: 2
    auto ofdmMod = new OFDMWithStaticPilot!(Complex!float)(8, 3, 4, 2);

    Complex!float[] inps = [Complex!float(1, 1), Complex!float(-1, -1), Complex!float(1, -1), Complex!float(-1, 1)];
    Complex!float[] res;

    // 変調
    ofdmMod.modulate(inps, res);

    assert(res.length == (8 + 3) * 2);

    // CPのチェック
    foreach(i; 0 .. 3 * 2)
        assert(res[i] == res[$ - 3*2 + i]);

    // 振幅変動と位相回転を与える
    foreach(ref e; res)
        e *= Complex!float(0.5, 1.0);

    Complex!float[] demodsym;
    ofdmMod.demodulate(res, demodsym);

    assert(demodsym.length == inps.length);

    foreach(i; 0 .. inps.length){
        assert(isClose(inps[i].re, demodsym[i].re));
        assert(isClose(inps[i].im, demodsym[i].im));
    }
}
