module dffdd.blockdiagram.quantizer;

import std.math;
import std.traits;
import std.range;
import std.complex;
import std.json;

import dffdd.utils.json;
import dffdd.utils.unit;


struct SimpleQuantizer(C)
{
    alias IType = C;
    alias OType = C;


    this(size_t nbit)
    {
        _nbit = nbit;
    }


    void opCall(IType input, ref OType output) @nogc
    {
        output = C(quantizeReal(input.re), quantizeReal(input.im));
    }


    void opCall(in IType[] input, OType[] output) @nogc
    in {
        assert(input.length == output.length);
    }
    do {
        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    typeof(this) dup() const pure nothrow @safe @nogc @property
    {
        return this;
    }


    JSONValue dumpInfoToJSON() const
    {
        return JSONValue([
            "bits": _nbit
        ]);
    }


  private:
    size_t _nbit;


    typeof(C.init.re) quantizeReal(typeof(C.init.re) x) @nogc
    {
        // 正・負の両方で_nbit-1ビットだけ精度を持つため，全体では_nbitの精度
        immutable long scale = 1L << (_nbit - 1);
        immutable real maxAmp = 30.dBm.volt;
        return floor(x / maxAmp * scale) / scale * maxAmp;
    }
}
