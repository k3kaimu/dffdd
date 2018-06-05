module dffdd.utils.unit;

import std.math;


// 系の基準となるインピーダンス
enum real IMPEDANCE = 1.0;


Gain dB(real db)
{
    return Gain(10.0^^(db/20));
}


struct Gain
{
  pure nothrow @safe @nogc
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
  }


  @property const pure nothrow @safe @nogc
  {
    deprecated
    real dB()
    {
        return 20 * log10(_g1V);
    }


    deprecated
    real gain() pure nothrow @safe @nogc
    {
        return _g1V;
    }


    real asV() pure nothrow @safe @nogc
    {
        return _g1V;
    }


    real asP() pure nothrow @safe @nogc
    {
       return _g1V^^2;
    }


    real asdB() pure nothrow @safe @nogc
    {
        return 20 * log10(_g1V);
    }
  }


  pure nothrow @safe @nogc 
  {
    Gain opBinary(string op : "*")(Gain g) const
    {
        return typeof(this)(_g1V * g._g1V);
    }


    Gain opBinary(string op : "/")(Gain g) const
    {
        return typeof(this)(_g1V / g._g1V);
    }


    void opOpAssign(string op)(Gain g)
    if(op == "*" || op == "/")
    {
        this = this.opBinary!op(g);
    }
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
    assert(approxEqual(20.dB.asV, 10));
    assert(approxEqual(Gain.fromPowerGain(100).asV, 10));
    assert(approxEqual(Gain.fromVoltageGain(123).asV, 123));
}



Voltage dBm(real dbm)
{
    return Voltage.fromdBm(dbm);
}


struct Voltage
{
    static
    Voltage fromdBm(real dbm)
    {
        return Voltage(sqrt(IMPEDANCE) * 10^^((dbm - 30)/20));
    }


    this(real v1V) pure nothrow @safe @nogc
    {
        _g1V = v1V;
    }


  @property const
  {
    real dBm() pure nothrow @safe @nogc
    {
        auto w = _g1V^^2 / IMPEDANCE;
        return 10*log10(w) + 30;
    }


    real volt() pure nothrow @safe @nogc
    {
        return _g1V;
    }


    real watt() pure nothrow @safe @nogc
    {
        return _g1V^^2 / IMPEDANCE;
    }
  }

    string toString() const
    {
        import std.format;
        return format("%sdBm", this.dBm);
    }


    Voltage opBinary(string op : "*")(Gain g) const pure nothrow @safe @nogc
    {
        return Voltage(this._g1V * g.asV);
    }


    Voltage opBinary(string op : "/")(Gain g) const pure nothrow @safe @nogc
    {
        return Voltage(this._g1V / g.asV);
    }


    Gain opBinary(string op : "/")(Voltage v) const pure nothrow @safe @nogc
    {
        return Gain.fromVoltageGain(this._g1V / v._g1V);
    }


    void opOpAssign(string op)(Gain g) pure nothrow @safe @nogc
    if(op == "*" || op == "/")
    {
        this = this.opBinary!op(g);
    }


  private:
    real _g1V;
}

unittest
{
    assert(approxEqual((-10).dBm.volt, sqrt(IMPEDANCE)*0.01));
    assert(approxEqual((10).dBm.volt, sqrt(IMPEDANCE)*0.1));
    assert(approxEqual((10).dBm.dBm, 10));
    assert(approxEqual((-10).dBm.dBm, -10));

    enum Voltage g1 = (-10).dBm;
    assert(approxEqual(g1.volt, sqrt(IMPEDANCE)*0.01));

    enum Voltage g2 = 10.dBm;
    assert(approxEqual(g2.volt, sqrt(IMPEDANCE)*0.1));
}
