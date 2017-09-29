module dffdd.blockdiagram.iqmixer;

import dffdd.utils.unit;
import dffdd.blockdiagram.noise;
import std.random;
import std.math;
import std.complex;
import std.range;


auto addIQImbalance(R)(R r, Gain gain, Gain irr)
{
    return IQImbalance!R(r, gain, irr);
}

struct IQImbalance(R)
{
    this(R r, Gain gain, Gain irr, real theta)
    {
        _r = r;
        _g1V = gain.gain;
        _g2V = gain.gain / irr.gain * std.complex.expi(theta);
    }


    auto front()
    {
        auto e = _r.front;
        return e * _g1V + e.conj * _g2V;
    }


    void popFront()
    {
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
    real _g1V;              // 電圧系での真値のゲイン
    Complex!real _g2V;      // 電圧計での真値のImageゲイン
}


struct IQImbalanceConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(Gain gain, Gain irr, real theta)
    {
        _g1V = gain.gain;
        _g2V = gain.gain / irr.gain * std.complex.expi(theta);
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        output = input * _g1V + input.conj * _g2V;
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        output.length = input.length;

        foreach(i; 0 .. input.length)
            this.opCall(input[i], output[i]);
    }


    typeof(this) dup() const pure nothrow @safe @nogc @property
    {
        return this;
    }


  private:
    real _g1V;      // 電圧系での真値のゲイン
    C _g2V;         // 電圧計での真値のImageゲイン
}

unittest
{
    import std.complex;
    import std.random;
    import std.algorithm;
    import dffdd.blockdiagram.utils;
    import std.math;

    alias C = Complex!float;

    auto signal = new C[64];
    foreach(ref e; signal)
        e = C(uniform01(), uniform01);

    auto r0 = IQImbalance!(C[])(signal, 10.dB, 10.dB, 0.1).array;
    auto r1 = signal.connectTo!(IQImbalanceConverter!C)(10.dB, 10.dB, 0.1);
    auto r2 = signal.chunks(3).connectTo!(IQImbalanceConverter!C)(10.dB, 10.dB, 0.1).joiner;

    assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r1));
    assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r2));
}


auto addPhaseNoise(R)(R r, Gain phaseNoise, real alpha = 0.5)
{
    return PhaseNoise!R(r, phaseNoise, alpha);
}


struct PhaseNoise(R)
{
    this(R r, real carrFreq, real samplingFreq, real paramC)
    {
        _r = r;
        _noiseGain = 4 * PI * carrFreq * sqrt(paramC / samplingFreq);
        _phi = 0;
        _noise = boxMullerNoise();
    }


    auto front()
    {
        return std.complex.expi(_phi) * _r.front;
    }


    void popFront()
    {
        _phi += _noise.front.re * _noiseGain;
        _noise.popFront();
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
        dst._noise = this._noise.save;

        return dst;
    }
  }


  private:
    R _r;
    real _noiseGain;
    real _phi;
    BoxMuller!Random _noise;
}




/+
struct PhaseNoise(C)
{
    this(Voltage phaseNoise, Gain iirAlpha)
    {
        _g1V = phaseNoise.V / SQRT2;
        _iir1 = IIRFilter1!C(iirAlpha.gain);
    }


    void apply(C1 : C)(in C1[] x, C1[] o)
    in{
        assert(x.length == o.length);
    }
    body{
        o[] = x[];
        
    }


  private:
    real _g1V;
    BoxMuler _noise;
    IIRFilter1!C _iir1;
}
+/