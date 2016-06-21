module dffdd.blockdiagram.iqmixer;

import dffdd.utils.unit;
import dffdd.blockdiagram.noise;
import std.random;
import std.math;


auto addIQImbalance(R)(R r, Gain gain, Gain irr)
{
    return IQImbalance!R(r, gain, irr);
}

struct IQImbalance(R)
{
    this(R r, Gain gain, Gain irr)
    {
        _r = r;
        _g1V = gain.gain;
        _g2V = gain.gain / irr.gain;
        //import std.stdio;
        //writefln("%s, %s", _g1V, _g2V);
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


  private:
    R _r;
    real _g1V;      // 電圧系での真値のゲイン
    real _g2V;      // 電圧計での真値のImageゲイン
}


auto addPhaseNoise(R)(R r, Gain phaseNoise, real alpha = 0.5)
{
    return PhaseNoise!R(r, phaseNoise, alpha);
}


struct PhaseNoise(R)
{
    this(R r, Gain phaseNoise, real alpha = 0.5)
    {
        _r = r;
        _g1V = phaseNoise.gain;
        _last = 0;
        _alpha = alpha;
        _noise = boxMullerNoise();
        _value = _noise.front.re * _g1V;
    }


    auto front()
    {
        return expi(_value) * _r.front;
    }


    void popFront()
    {
        _last = _value;
        _noise.popFront();
        _r.popFront();

        auto now = _noise.front.re * _g1V;
        _value = now - _last * _alpha;
    }


    bool empty()
    {
        return _r.empty;
    }


  private:
    R _r;
    real _g1V;
    real _last;
    real _value;
    real _alpha;
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