module dffdd.mod.qpsk;

import std.bitmanip;
import std.complex;


struct QPSK
{
    alias InputElementType = ubyte;
    alias OutputElementType = Complex!float;

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
            outputs[i] = Complex!float(inputs[i*2+0] != 0 ? -SQRT1_2 : SQRT1_2,
                                       inputs[i*2+1] != 0 ? -SQRT1_2 : SQRT1_2);

        return outputs;
    }


    static
    ref OutputElementType[] modulate(const ref BitArray inputs, return ref OutputElementType[] outputs)
    in{
        assert(inputs.length % symInputLength == 0);
    }
    body{
        import std.math;

        if(outputs.length != inputs.length/2)
            outputs.length = inputs.length/2;

        foreach(i; 0 .. inputs.length/2)
            outputs[i] = Complex!float(inputs[i*2+0] ? -SQRT1_2 : SQRT1_2,
                                       inputs[i*2+1] ? -SQRT1_2 : SQRT1_2);

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
    ref BitArray demodulate(in OutputElementType[] inputs, return ref BitArray outputs)
    {
        if(outputs.length != inputs.length*2)
            outputs.length = inputs.length*2;

        foreach(i; 0 .. inputs.length){
            outputs[i*2+0] = inputs[i].re < 0;
            outputs[i*2+1] = inputs[i].im < 0;
        }

        return outputs;
    }
}

unittest
{
    BitArray bits;
    foreach(i; 0 .. 52)
        bits ~= i % 3 == 0;

    Complex!float[] signal;
    BitArray copy = bits;
    QPSK qpsk;
    qpsk.modulate(bits, signal);
    assert(signal.length == 26);
    bits.length = 0;
    qpsk.demodulate(signal, bits);

    assert(bits.length == 52);
    foreach(i; 0 .. 52)
        assert(bits[i] == copy[i]);
}


real berQPSKFromSNR(real snr)
{
    import dffdd.mod.qam;

    return berQAMFromSNR(snr, 4);
}
