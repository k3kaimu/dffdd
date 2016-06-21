module dffdd.mod.ofdm;


import std.stdio;

import dffdd.itpp;


struct OFDM
{
    alias InputElementType = cfloat;
    alias OutputElementType = InputElementType;


    this(uint nFFT, uint nCp, uint nTone, uint nUpSampling = 1)
    {
        _impl = makeOFDMModulator(nFFT, nCp, nUpSampling);
        _nFFT = nFFT;
        _nCp = nCp;
        _nTone = nTone;
        _nUpSampling = nUpSampling;
    }


    size_t symInputLength() const @property { return _nTone; }
    size_t symOutputLength() const @property { return (_nFFT + _nCp) * _nUpSampling; }


    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    in{
        assert(inputs.length % this.symInputLength == 0);
        assert(outputs.length % this.symOutputLength == 0);
    }
    body{
        InputElementType[] buf = new InputElementType[_nFFT];

        outputs.length = inputs.length / _nTone * this.symOutputLength;
        foreach(i; 0 .. inputs.length / _nTone){
            buf[] = 0+0i;
            if(_nTone == _nFFT){
                buf[] = inputs[i*_nTone .. (i+1)*_nTone];
            }else if(_nTone == _nFFT -1){
                buf[1 .. $] = inputs[i*_nTone .. (i+1)*_nTone];
            }else if(_nTone % 2 == 0){
                buf[1 .. _nTone/2+1] = inputs[i*_nTone .. i*_nTone + _nTone/2];
                buf[$ - _nTone/2 .. $] = inputs[i*_nTone + _nTone/2 .. (i+1)*_nTone];
            }else{
                buf[1 .. (_nTone+1)/2+1] = inputs[i*_nTone .. i*_nTone + (_nTone+1)/2];
                buf[$ - _nTone/2 .. $] = inputs[i*_nTone + (_nTone+1)/2 .. (i+1)*_nTone];
            }

            auto dst = outputs[i*symOutputLength .. (i+1)*symOutputLength];
            _impl.modulate(buf, dst);
            assert(dst.ptr == outputs.ptr + i*symOutputLength);
        }

        return outputs;
    }


    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    in{
        assert(inputs.length % this.symOutputLength == 0);
        assert(outputs.length % this.symInputLength == 0);
    }
    body{
        InputElementType[] buf = new InputElementType[_nFFT];
        immutable slen = this.symOutputLength;

        outputs.length = this.symInputLength * inputs.length / this.symOutputLength;

        foreach(i; 0 .. inputs.length / slen) {
            _impl.demodulate(inputs[i*slen .. (i+1)*slen], buf);

            if(_nTone == _nFFT){
                outputs[i*_nTone .. (i+1)*_nTone] = buf[];
            }else if(_nTone == _nFFT - 1){
                outputs[i*_nTone .. (i+1)*_nTone] = buf[1 .. $];
            }else if(_nTone % 2 == 0){
                outputs[i*_nTone .. i*_nTone + _nTone/2] = buf[1 .. _nTone/2+1];
                outputs[i*_nTone + _nTone/2 .. (i+1)*_nTone] = buf[$ - _nTone/2 .. $];
            }else{
                outputs[i*_nTone .. i*_nTone + (_nTone+1)/2] = buf[1 .. (_nTone+1)/2+1];
                outputs[i*_nTone + (_nTone+1)/2 .. (i+1)*_nTone] = buf[$ - _nTone/2 .. $];
            }
        }

        return outputs;
    }


  private:
    typeof(makeOFDMModulator(1, 1, 1)) _impl;
    uint _nFFT;
    uint _nCp;
    uint _nTone;
    uint _nUpSampling;
}

