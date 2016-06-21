module dffdd.utils.unit;

import std.math;


Gain dB(real db)
{
    return Gain(10.0^^(db/20));
}


Gain gain(real g1V)
{
    return Gain(g1V);
}


struct Gain
{
    this(real g1V)
    {
        _g1P = sqrt(g1V);
    }


  @property
  {
    real dB() pure nothrow @safe @nogc
    {
        return 10 * log10(_g1P);
    }


    real gain() pure nothrow @safe @nogc
    {
        return _g1P ^^ 2;
    }


    real gainP() pure nothrow @safe @nogc
    {
        return _g1P;
    }
  }


  private:
    real _g1P;
}



Voltage dBm(real dbm)
{
    return Voltage(10.0^^((dbm - 30)/20));
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


  @property
  {
    real dBm() pure nothrow @safe @nogc
    {
        // 1V => 1W == 30dBm
        return 30 + 20*log10(_g1V);
    }


    real dBW() pure nothrow @safe @nogc
    {
        return 20*log10(_g1V);
    }


    real V() pure nothrow @safe @nogc
    {
        return _g1V;
    }


    real W() pure nothrow @safe @nogc
    {
        return _g1V^^2;
    }
  }


  private:
    real _g1V;
}