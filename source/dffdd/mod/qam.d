module dffdd.mod.qam;

import dffdd.itpp;
import std.complex;
import std.mathspecial;
import std.math;
import std.numeric;


struct QAM
{
    alias InputElementType = ubyte;
    alias OutputElementType = Complex!float;


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


real ber16QAMFromSNR(real snr)
{
    return 4.0L / 4.0L * (1.0L - 1/4.0L) / 2 * erfc(sqrt(3.0L * 4 / 2 / (16 - 1) * 10.0L^^((snr-6)/10)));
}


real snrFromBer16QAM(real ber)
{
    return findRoot(delegate(real snr){ return (ber16QAMFromSNR(snr) - ber) / ber; }, 0.0L, 25.0L);
}
