module dffdd.mod.qam;

import dffdd.itpp;

struct QAM
{
    alias InputElementType = ubyte;
    alias OutputElementType = cfloat;


    this(uint mary)
    {
        _impl = makeQAMModulator(mary);
        _mary = mary;
        while(mary != 1){
            mary >>= 1;
            ++_nIL;
        }
    }


    size_t symInputLength() const @property { return _nIL; }
    size_t symOutputLength() const @property { return 1; }


    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    {
        outputs = _impl.modulate_bits(inputs, outputs);
        return outputs;
    }


    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    {
        outputs = _impl.demodulate_bits(inputs, outputs);
        return outputs;
    }


  private:
    typeof(makeQAMModulator(1)) _impl;
    uint _mary;
    uint _nIL;
}