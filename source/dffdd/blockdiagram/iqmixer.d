module dffdd.blockdiagram.iqmixer;

import dffdd.utils.unit;
import dffdd.blockdiagram.noise;
import dffdd.utils.json;
import std.random;
import std.math;
import std.complex;
import std.range;
import std.traits;
import std.json;


auto addIQImbalance(R)(R r, Gain gain, Gain irr)
{
    return IQImbalance!R(r, gain, irr);
}


struct IQImbalanceConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;
    alias F = typeof(C.init.re);


    this(Gain gain, Gain irr, real theta)
    {
        this(gain, cast(C)(1.0L / irr.asV  * std.complex.expi(theta)));
    }


    this(Gain gain, C coef)
    {
        _g1V = gain.asV;
        _g2V = gain.asV * coef;
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        output = input * _g1V + input.conj * _g2V;
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        if(output.length != input.length)
            output.length = input.length;

        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    C imbCoef() const @property
    {
        return _g2V / _g1V;
    }


    Gain gain() const @property
    {
        return Gain.fromPowerGain(_g1V^^2 + _g2V.sqAbs());
    }


    JSONValue dumpInfoToJSON() const
    {
        JSONValue jv = JSONValue(string[string].init);
        jv["gain"] = DefaultJSONEnv.toJSONValue(_g1V);
        jv["imbCoef"] = DefaultJSONEnv.toJSONValue(_g2V / _g1V);
        return jv;
    }


    typeof(this) dup() const pure nothrow @safe @nogc @property
    {
        return this;
    }


  private:
    F _g1V;      // 電圧系での真値のゲイン
    C _g2V;         // 電圧計での真値のImageゲイン
}


struct IQImbalance
{
    static
    auto makeBlock(R)(R range, Gain gain, Complex!real coef)
    {
        import dffdd.blockdiagram.utils : connectTo;

        alias E = Unqual!(ElementType!R);
        return range.connectTo!(IQImbalanceConverter!E)(gain, coef);
    }
}


unittest 
{
    import std.algorithm : map, equal;
    import dffdd.blockdiagram.utils;

    static assert(isDuplicatableConverter!(IQImbalanceConverter!(Complex!real)));

    Complex!real[] signal = new Complex!real[1024];
    foreach(i; 0 .. 1024)
        signal[i] = Complex!real(uniform01, uniform01);
    enum coef = Complex!real(1.0, 2.0);

    auto iqimb = signal.connectTo!IQImbalance(0.dB, coef);
    assert(equal(iqimb, signal.map!(a => a + a.conj * Complex!real(1.0, 2.0))));
}


auto addPhaseNoise(R)(R r, Gain phaseNoise, real alpha = 0.5)
{
    return PhaseNoise!R(r, phaseNoise, alpha);
}


deprecated
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


struct FROPhaseNoiseGenerator(C, Rnd)
{
    this(Rnd rndGen, real sampFreqHz, real betaBWHz)
    {
        _normGen = BoxMuller!Rnd(rndGen);
        _phase = 0;
        _scale = 2 * sqrt(PI * betaBWHz / sampFreqHz);
    }


    void opCall(in C[] input, ref C[] output)
    {
        if(output.length != input.length)
            output.length = input.length;
        
        foreach(i; 0 .. input.length) {
            output[i] = input[i] * std.complex.expi(_phase);
            _phase += _normGen.front.re * _scale;
            _normGen.popFront();
        }
    }


    FROPhaseNoiseGenerator!(C, Rnd) dup()
    {
        FROPhaseNoiseGenerator!(C, Rnd) dst;
        dst._phase = this._phase;
        dst._scale = this._scale;
        dst._normGen = this._normGen.save;
        return dst;
    }


  private:
    alias F = typeof(C.init.re);

    F _phase;
    F _scale;
    BoxMuller!Rnd _normGen;
}



// struct PLLPhaseNoiseGenerator(C, Rnd)
// {
//     this(Rnd rnd, real Cvco, real Cxtl, real kdp, real omegaLp, real sampFreqHz, real ccontr, real fc)
// }




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