module dffdd.mod.qpsk;

import std.complex;

import dffdd.mod.primitives : Bit;

struct QPSK(C)
{
    alias InputElementType = Bit;
    alias OutputElementType = C;

    enum size_t symInputLength = 2;
    enum size_t symOutputLength = 1;


    static
    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    in{
        assert(inputs.length % symInputLength == 0);
    }
    body{
        import std.math;

        if(outputs.length != inputs.length/2)
            outputs.length = inputs.length/2;

        foreach(i; 0 .. inputs.length/2)
            outputs[i] = C(inputs[i*2+0] != 0 ? -SQRT1_2 : SQRT1_2,
                           inputs[i*2+1] != 0 ? -SQRT1_2 : SQRT1_2);

        return outputs;
    }


    static
    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    {
        if(outputs.length != inputs.length*2)
            outputs.length = inputs.length*2;

        foreach(i; 0 .. inputs.length){
            outputs[i*2+0] = inputs[i].re > 0 ? 0 : 1;
            outputs[i*2+1] = inputs[i].im > 0 ? 0 : 1;
        }

        return outputs;
    }


    static
    ref ushort[] demodulate_symbol(in OutputElementType[] inputs, return ref ushort[] symbols)
    {
        if(symbols.length != inputs.length)
            symbols.length = inputs.length;

        foreach(i; 0 .. inputs.length) {
            symbols[i] = 0;
            symbols[i] += inputs[i].re > 0 ? 0 : 2;
            symbols[i] += inputs[i].im > 0 ? 0 : 1;
        }

        return symbols;
    }
}


real berQPSKFromSNR(real snr)
{
    import dffdd.mod.qam;

    return berQAMFromSNR(snr, 4);
}
