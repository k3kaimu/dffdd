module dffdd.itpp;


import std.stdio;
import std.algorithm;
import std.typecons;
import std.range;
import std.random;


extern(C++, itpp)
{
    interface QAM{}
    interface OFDM{}
}


extern(C++, WITPP)
{
    extern(C++, QAM)
    {
        itpp.QAM newObject(uint mary);
        void deleteObject(itpp.QAM);

        uint demodulate_bits(const itpp.QAM, const(float)*, uint, ubyte*, uint*);
        uint modulate_bits(const itpp.QAM, const(ubyte)*, uint, float*, uint*);
    }

    extern(C++, OFDM)
    {
        itpp.OFDM newObject(uint inNfft, uint inNcp, uint inNupsample);
        void deleteObject(itpp.OFDM);

        uint demodulate(itpp.OFDM, const(float)*, uint, float*, uint*);
        uint modulate(itpp.OFDM, const(float)*, uint, float*, uint*);
    }
}


auto makeQAMModulator(uint mary)
{
    static struct QAMModulator
    {
        this(itpp.QAM qam) { _qam = qam; }

        ~this()
        {
            if(_qam !is null)
                WITPP.QAM.deleteObject(_qam);

            _qam = null;
        }


        ref cfloat[] modulate_bits(in ubyte[] bits, ref return cfloat[] outputs) const
        {
            auto size = WITPP.QAM.modulate_bits(_qam, cast(const(ubyte)*)bits.ptr, cast(uint)bits.length, null, null);
            outputs.length = size;
            WITPP.QAM.modulate_bits(null, null, 0, cast(float*)outputs.ptr, &size);
            assert(size == outputs.length);

            return outputs;
        }


        ref ubyte[] demodulate_bits(in cfloat[] inputs, ref return ubyte[] bits) const
        {
            auto size = WITPP.QAM.demodulate_bits(_qam, cast(const(float)*)inputs.ptr, cast(uint)inputs.length, null, null);
            bits.length = size;
            WITPP.QAM.demodulate_bits(null, null, 0, bits.ptr, &size);
            assert(size == bits.length);

            return bits;
        }


      private:
        itpp.QAM _qam;
    }


    return RefCounted!QAMModulator(WITPP.QAM.newObject(mary));
}


auto makeOFDMModulator(uint inNfft, uint inNcp, uint inNupsample)
{
    static struct OFDMModulator
    {
        this(itpp.OFDM ofdm) { _ofdm = ofdm; }

        ~this()
        {
            if(_ofdm !is null)
                WITPP.OFDM.deleteObject(_ofdm);

            _ofdm = null;
        }


        ref cfloat[] modulate(in cfloat[] inputs, ref return cfloat[] outputs)
        {
            auto size = WITPP.OFDM.modulate(_ofdm, cast(const(float)*)inputs.ptr, cast(uint)inputs.length, null, null);
            outputs.length = size;
            WITPP.OFDM.modulate(null, null, 0, cast(float*)outputs.ptr, &size);
            assert(size == outputs.length);

            return outputs;
        }


        ref cfloat[] demodulate(in cfloat[] inputs, ref return cfloat[] outputs)
        {
            auto size = WITPP.OFDM.demodulate(_ofdm, cast(const(float)*)inputs.ptr, cast(uint)inputs.length, null, null);
            outputs.length = size;
            WITPP.OFDM.demodulate(null, null, 0, cast(float*)outputs.ptr, &size);
            assert(size == outputs.length);

            return outputs;
        }


      private:
        itpp.OFDM _ofdm;
    }


    return RefCounted!OFDMModulator(WITPP.OFDM.newObject(inNfft, inNcp, inNupsample));
}
