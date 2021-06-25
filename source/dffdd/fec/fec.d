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
