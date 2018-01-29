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
        this(r, gain, (1.0L/irr.gain) * std.complex.expi(theta));
    }


    this(R r, Gain gain, Complex!real imbCoef)
    {
        _r = r;
        _g1V = gain.gain;
        _g2V = gain.gain * imbCoef;
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
    real _g1V;      // 電圧系での真値のゲイン
    Complex!real _g2V;      // 電圧計での真値のImageゲイン
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