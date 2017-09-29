module dffdd.utils.unit;

import std.math;


private
bool isNaNCTFE(real x)
{
    return x != x;
}

unittest
{
    assert(!isNaNCTFE(1.0L));
    assert(!isNaNCTFE(real.infinity));
    assert(isNaNCTFE(real.nan));
}


private
bool isInfinityCTFE(real x)
{
    return !isNaNCTFE(x) && x - x != 0;
}

unittest
{
    assert(!isInfinityCTFE(1.0L));
    assert(!isInfinityCTFE(real.nan));
    assert(isInfinityCTFE(real.infinity));
    assert(isInfinityCTFE(-real.infinity));
}


private
long floorStupid(real x)
{
    if(x % 1 == 0)
        return cast(long)x;
    else{
        return (cast(long)x) + (x < 0 ? -1 : 0);
    }
}

unittest
{
    assert(floorStupid(0) == 0);
    assert(floorStupid(0.999) == 0);
    assert(floorStupid(1) == 1);
    assert(floorStupid(-1) == -1);
    assert(floorStupid(-1.00001) == -2);
    assert(floorStupid(-2.99999) == -3);
}


private
real expCTFE(real x)
{
    if(isNaNCTFE(x))
        return real.nan;
    else if(x < 0)
        return 1.0L / (expCTFE(-x));
    else if(isInfinityCTFE(x) || x > int.max)
        return real.infinity;
    else if(x <= 1){
        real r = 0;
        real d = 1;

        // マクローリン展開
        foreach(i; 0 .. 10){
            if(i != 0){
                d *= x;
                d /= i;
            }
            r += d;
        }

        return r;
    }else{
        immutable long n = floorStupid(x);
        real r = 1;

        // 整数部分の計算
        foreach(i; 0 .. n)
            r *= E;

        // 整数部分と小数点以下の結果を掛け合わせる
        return r * expCTFE(x - n);
    }
}

unittest
{
    foreach(i; 0 .. 2000){
        real r = i / 10.0L - 100;
        assert(approxEqual(expCTFE(r), exp(r)));
    }
}


Gain dB(real db)
{
    if(__ctfe)
        return Gain(expCTFE(db / 20 * log(10.0L)));
    else
        return Gain(10.0^^(db/20));
}


struct Gain
{
    static
    Gain fromPowerGain(real g)
    {
        return Gain(sqrt(g));
    }


    static
    Gain fromVoltageGain(real g)
    {
        return Gain(g);
    }


    private this(real g1V)
    {
        _g1V = g1V;
    }


  @property const
  {
    real dB() pure nothrow @safe @nogc
    {
        return 20 * log10(_g1V);
    }


    deprecated
    real gain() pure nothrow @safe @nogc
    {
        return _g1V;
    }


    // real gainP() pure nothrow @safe @nogc
    // {
    //    return _g1P;
    // }
  }


    string toString() const
    {
        import std.format;
        return format("%sdB", log10(_g1V)*20);
    }


  private:
    real _g1V;
}

unittest
{
    import std.stdio;
    assert(approxEqual(20.dB.gain, 10));
    assert(approxEqual(Gain.fromPowerGain(100).gain, 10));
    assert(approxEqual(Gain.fromVoltageGain(123).gain, 123));
}



Voltage dBm(real dbm)
{
    if(__ctfe)
        return Voltage(sqrt(50.0) * expCTFE((dbm - 30)/20 * log(10.0)));
    else
        return Voltage(sqrt(50.0) * 10^^((dbm - 30)/20));
}


Voltage V(real v)
{
    return Voltage(v);
}



struct Voltage
{
    this(real v1V)
    {
        _g1V = v1V;
    }


  @property const
  {
    real dBm() pure nothrow @safe @nogc
    {
        auto w = _g1V^^2 / 50;
        return 10*log10(w) + 30;
    }


    //deprecated
    //real dBW() pure nothrow @safe @nogc
    //{
    //    return 20*log10(_g1V);
    //}


    deprecated
    real V() pure nothrow @safe @nogc
    {
        return _g1V;
    }


    //deprecated
    //real W() pure nothrow @safe @nogc
    //{
    //    return _g1V^^2;
    //}
  }

    string toString() const
    {
        import std.format;
        return format("%sdBm", this.dBm);
    }


  private:
    real _g1V;
}

unittest
{
    assert(approxEqual((-10).dBm.V, sqrt(50.0)*0.01));
    assert(approxEqual((10).dBm.V, sqrt(50.0)*0.1));
    assert(approxEqual((10).dBm.dBm, 10));
    assert(approxEqual((-10).dBm.dBm, -10));

    enum Voltage g1 = (-10).dBm;
    assert(approxEqual(g1.V, sqrt(50.0)*0.01));

    enum Voltage g2 = 10.dBm;
    assert(approxEqual(g2.V, sqrt(50.0)*0.1));
}
