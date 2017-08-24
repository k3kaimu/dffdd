module dffdd.mod.sccp;


import dffdd.utils.fft;


final class SCCP(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;
    private alias F = typeof(C.init.re);


    this(uint nBlock, uint nCp, uint nOS)
    {
        auto resp = new C[nBlock * nOS];
        foreach(ref e; resp)
            e = C(1);
        
        this(nBlock, nCp, nOS, resp);
    }


    this(uint nBlock, uint nCp, uint nOS, const(C)[] freqResp)
    {
        _nBlock = _nBlock;
        _nCp = nCp;
        _nOS = nOS;
        _freqResp = freqResp.dup;
    }


    size_t symInputLength() const @property { return _nBlock; }
    size_t symOutputLength() const @property { return (_nBlock + _nCp) * _nOS; }

    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    in{
        assert(inputs.length % this.symInputLength == 0);
    }
    body{
        if(outputs.length != inputs.length / this.symInputLength * this.symOutputLength)
            outputs.length = inputs.length / this.symInputLength * this.symOutputLength;

        auto cp = outputs[0 .. _nCP * _nOS];
        auto sy = outputs[_nCP * _nOS];

        foreach(i; 0 .. inputs.length){
            foreach(r; 0 .. _nOS)
                outputs[i*_nOS + r] = inputs[i];
        }

        // add CP
        cp[] = sy[$ - cp.length .. $];

        return outputs;
    }


    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    in{
        assert(inputs.length % this.symOutputLength == 0);
    }
    body{
        if(outputs.length != inputs.length / this.symOutputLength * this.symInputLength)
            outputs.length = inputs.length / this.symOutputLength * this.symInputLength;

        auto fftwInputs = _fftw.inputs!F;
        auto fftwOutputs = _fftw.outputs!F;

        foreach(i; 0 .. inputs.length / this.symOutputLength)
        {
            // remove CP
            auto sy = inputs[i * this.symOutputLength .. (i+1) * this.symOutputLength][_nCP * _nOS .. $];

            fftwInputs[] = sy[];
            _fftw.fft!F();

            foreach(i; 0 .. _nBlock * _nOS)
                fftwInputs[i] = fftwOutputs[i] * _freqResp[i];

            _fftw.ifft!F();
            outputs[i * this.symInputLength .. (i+1) * this.symInputLength] = fftwOutputs[];
        }

        return outputs;
    }


    inout(C)[] frequencyResponse() inout @property
    {
        return _freqResp;
    }


  private:
    uint _nBlock;
    uint _nCp;
    uint _nOS;
    C[] _freqResp;
}
