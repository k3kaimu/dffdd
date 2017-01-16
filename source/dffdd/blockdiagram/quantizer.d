module dffdd.blockdiagram.quantizer;

import std.math;
import std.traits;
import std.range;
import std.complex;


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


struct SimpleQuantizer(R)
{
    alias E = typeof(Unqual!(ElementType!R).init.re);

    this(R r, size_t nbit)
    {
        _r = r;
        _nbit = nbit;
    }


    Complex!E front() @property
    {
        auto f = _r.front * (1 << (_nbit - 1));
        long ivr, ivi;

        if(f.re >= long.max)
            ivr = long.max;
        else if(f.re <= long.min)
            ivr = long.min;
        else
            ivr = cast(long)f.re;

        if(f.im >= long.max)
            ivi = long.max;
        else if(f.im <= long.min)
            ivi = long.min;
        else
            ivi = cast(long)f.im;

        return Complex!E(ivr, ivi) / (1 << (_nbit - 1));
    }


    void popFront()
    {
        _r.popFront();
    }


  static if(isInfinite!R)
    enum bool empty = true;
  else
  {
    bool empty() @property
    {
        return _r.empty;
    }
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
}
