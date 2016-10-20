module dffdd.mod.bpsk;

import std.complex;
import std.mathspecial;
import std.math;
import std.numeric;


struct BPSK
{
    alias InputElementType = ubyte;
    alias OutputElementType = Complex!float;

    enum size_t symInputLength = 1;
    enum size_t symOutputLength = 1;


    static
    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    {
        if(outputs.length != inputs.length)
            outputs.length = inputs.length;

        foreach(i; 0 .. inputs.length)
            outputs[i] = Complex!float(inputs[i] != 0 ? -1 : 1, 0);

        return outputs;
    }


    static
    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    {
        if(outputs.length != inputs.length)
            outputs.length = inputs.length;

        foreach(i; 0 .. inputs.length)
            outputs[i] = inputs[i].re > 0 ? 0 : 1;

        return outputs;
    }
}


real berBPSKFromSNR(real snr)
{
    import dffdd.mod.qam;

    return berQAMFromSNR(snr + 10*log10(2.0), 4);
}
