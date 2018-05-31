module dffdd.blockdiagram.quantizer;

import std.math;
import std.traits;
import std.range;
import std.complex;
import std.json;

import dffdd.utils.json;
import dffdd.utils.unit;


struct Quantizer(R)
{
    import std.stdio;

    this(R r, size_t nbit)
    {
        _r = r;
        _nbit = nbit;
        _avgP = 0;
        _gain1V = 1;
        _i = 0;
    }


    auto front()
    {
        auto e = _r.front * _gain1V / 2 + 0.5+0.5i;
        auto re = e.re * (2.0^^_nbit - 1);
        auto im = e.im * (2.0^^_nbit - 1);
        auto sre = cast(long)re;
        auto sim = cast(long)im;

        auto ret = ((sre / (2.0^^_nbit - 1) + sim / (2.0^^_nbit - 1) * 1i) - (0.5+0.5i))*2 / _gain1V;
        return ret;
    }


    void popFront()
    {
        auto e = _r.front;
        auto p = e.re^^2 + e.im^^2;
        _avgP += p;
        ++_i;

        if(_i == 1024){
            _gain1V = _avgP == 0 ? 0 : 1/sqrt(_avgP / 1024) / sqrt(10.0f);     // PAPR10dBを想定
            _i = 0;
            _avgP = 0;
        }

        _r.popFront();
    }


    bool empty()
    {
        return _r.empty;
    }


  static if(isForwardRange!R)
  {
    typeof(this) save() @property
    {
        typeof(return) dst = this;

        dst._r = this._r.save;

        return dst;
    }
  }



  private:
    R _r;
    size_t _nbit;
    real _avgP;
    real _gain1V;
    size_t _i;
}


struct SimpleQuantizerConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(size_t nbit)
    {
        _nbit = nbit;
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        output = C(quantizeReal(input.re), quantizeReal(input.im));
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        output.length = input.length;

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


    typeof(C.init.re) quantizeReal(typeof(C.init.re) x)
    {
        // 正・負の両方で_nbit-1ビットだけ精度を持つため，全体では_nbitの精度
        immutable long scale = 1L << (_nbit - 1);
        immutable real maxAmp = 30.dBm.volt;
        return floor(x / maxAmp * scale) / scale * maxAmp;
    }
}


struct SimpleQuantizer
{
    static
    auto makeBlock(R)(R r, size_t nbit)
    {
        import dffdd.blockdiagram.utils : connectTo;
        alias E = Unqual!(ElementType!R);
        return r.connectTo!(SimpleQuantizerConverter!E)(nbit);
    }
}

unittest 
{
    auto r = SimpleQuantizer.makeBlock(Complex!float[].init, 1);
}
