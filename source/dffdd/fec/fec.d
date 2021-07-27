module dffdd.fec.fec;

import std.math;
import std.complex;

import dffdd.mod.primitives : Bit;

interface Fec(T)
{
    Bit[] encode(in Bit[] info, return ref Bit[] encoded);
    Bit[] decode(in T[] input, return ref Bit[] decoded);
}


interface BlockFec(T) : Fec!(T)
{
    size_t inputLength() const;
    size_t outputLength() const;


    override
    Bit[] encode(in Bit[] info, return ref Bit[] encoded)
    in {
        assert(info.length % this.inputLength == 0);
    }
    out(r) {
        // assert(r.length % this.outputLength == 0);
        // assert(info.length / this.inputLength == r.length / this.outputLength);
    };


    override
    Bit[] decode(in T[] input, return ref Bit[] decoded)
    in {
        assert(input.length % this.outputLength == 0);
    }
    out(r) {
        // assert(r.length % this.inputLength == 0);
        // assert(input.length / this.outputLength == r.length / this.inputLength);
    };
}



class NopBlockFec(T, alias decision) : Fec!T
if(is(typeof(decision(T.init) ) == Bit))
{
    this(size_t N)
    {
        _size = N;
    }


    size_t inputLength() const { return _size; }
    size_t outputLength() const { return _size; }


    Bit[] encode(in Bit[] info, return ref Bit[] encoded)
    {
        encoded.length = info.length;
        foreach(i; 0 .. info.length)
            encoded[i] = info[i];

        return encoded;
    }


    Bit[] decode(in T[] input, return ref Bit[] decoded)
    {
        decoded.length = input.length;
        foreach(i; 0 .. input.length)
            decoded[i] = decision(input[i]);
        
        return decoded;
    }


  private:
    size_t _size;
}
