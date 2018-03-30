module dffdd.blockdiagram.amplifier;

import std.traits;
import std.range;
import std.math;
import std.complex;
import std.json;

import dffdd.utils.unit;
import dffdd.utils.json;


struct PowerAmplifier(R)
{
    this(R r, Gain gain, Voltage iip3, Voltage iip5 = Voltage(0)/*, Voltage iip7 = Voltage(0)*/)
    {
        _r = r;

        _g1V = gain.gain;
        _g3V = iip3.V == 0 ? 0 : (gain.gain / iip3.V^^2);
        _g5V = iip5.V == 0 ? 0 : (gain.gain / iip5.V^^4);
        //_g7V = iip7.V == 0 ? 0 : (gain.gain / iip7.V^^6);
    }


    auto front()
    {
        auto x = _r.front;
        auto x1p = x.re^^2 + x.im^^2,
             x3 = x * x1p,
             x5 = x3 * x1p;
        //     x7 = x5 * x1p;

        return x * _g1V + x3 * _g3V + x5 * _g5V;// + x7 * _g7V;
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
        typeof(this) dst;
        dst._r = this._r.save;
        dst._g1V = this._g1V;
        dst._g3V = this._g3V;
        dst._g5V = this._g5V;

        return dst;
    }
  }


  private:
    R _r;
    real _g1V;
    real _g3V;
    real _g5V;
    //real _g7V;
}


/**
o: saturation value
s: smooth factor
u: |x[n]|
g(u) = (G*u)/(1+(u/o)^^(2s))^^(1/(2*s))
*/
struct RappModelConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;

    this(Gain gain, real smoothFactor, real saturation)
    {
        _g = gain.gain;
        _s = smoothFactor;
        _o = saturation;
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        auto r = std.complex.abs(input),
             u = input / r;         // unit vector

        // rが小さすぎるときに，単位ベクトルが発散するのを防ぐ
        if(r <= 1E-6){
            output = input;
        }else{
            r = (r * _g) / (( 1 + (r/_o)^^(2*_s) )^^(1/(2*_s)));
            output = r * u;
        }
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        output.length = input.length;

        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    /**
    入力信号が複素正規分布に従うと仮定して，出力電力を算出する
    */
    Voltage outputVoltage(Voltage x) const
    {
        return x;
    }


    typeof(this) dup() const pure nothrow @safe @nogc @property
    {
        return this;
    }


    JSONValue dumpInfoToJSON() const
    {
        return JSONValue([
            "gain":         _g,
            "smoothness":   _s,
            "saturation":   _o
        ]);
    }


  private:
    real _g, _s, _o;
}


struct RappModel
{
    static
    auto makeBlock(R)(R range, Gain gain, real smoothFactor, real saturation)
    {
        import dffdd.blockdiagram.utils : connectTo;
        alias E = Unqual!(ElementType!R);
        return range.connectTo!(RappModelConverter!E)(gain, smoothFactor, saturation);
    }
}

unittest
{
    auto r = RappModel.makeBlock(Complex!float[].init, 1.0.dB, 1, 1);
}


struct VGAConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;

    this(Gain gain)
    {
        _gain1V = gain.gain;
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        output = input * _gain1V;
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        output.length = input.length;

        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    Gain gain() const @property
    {
        return Gain.fromVoltageGain(_gain1V);
    }


    JSONValue dumpInfoToJSON() const
    {
        JSONValue jv = JSONValue(string[string].init);
        jv["gain"] = DefaultJSONEnv.toJSONValue(_gain1V);
        return jv;
    }


    VGAConverter!C dup() const pure nothrow @safe @nogc @property 
    {
        return this;
    }


  private:
    real _gain1V;
}


struct VGA(R)
{
    this(R r, Gain gain)
    {
        _r = r;
        _gain1V = gain.gain;
    }


    auto front()
    {
        return _r.front * _gain1V;
    }


    bool empty()
    {
        return _r.empty;
    }


    void popFront()
    {
        _r.popFront();
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
    real _gain1V;
}


struct PowerControlAmplifierConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(Voltage op, size_t avgSize = 512, size_t totalCount = 30)
    {
        _power = op.volt^^2;
        _alpha = 1;
        _cnt = 0;
        _avgSize = avgSize;
        _sumPower = 0;
        _avgCount = 0;
        _totalCount = totalCount;
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        _sumPower += input.re^^2 + input.im^^2;
        ++_cnt;
        if(_cnt == _avgSize){
            if(_avgCount < _totalCount){
                if(_sumPower == 0)
                    _alpha = _alpha;
                else
                    _alpha = sqrt(_power / (_sumPower / _avgSize));

                ++_avgCount;
            }

            _sumPower = 0;
            _cnt = 0;
        }

        output = input * _alpha;
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        output.length = input.length;

        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    Gain gain() @property const
    {
        return Gain.fromVoltageGain(_alpha);
    }


    PowerControlAmplifierConverter!C dup() const pure nothrow @safe @nogc @property 
    {
        return this;
    }


    JSONValue dumpInfoToJSON() const
    {
        return JSONValue([
            "gain": _alpha
        ]);
    }


  private:
    real _power;
    real _alpha;
    size_t _cnt;
    size_t _avgSize;
    size_t _avgCount;
    size_t _totalCount;
    real _sumPower;
}



struct PowerControlAmplifier
{
    static
    auto makeBlock(R)(R r, Voltage op, size_t avgSize = 512, size_t totalCount = 30)
    {
        import dffdd.blockdiagram.utils : connectTo;
        alias E = Unqual!(ElementType!R);
        return r.connectTo!(PowerControlAmplifierConverter!E)(op, avgSize, totalCount);
    }
}

unittest
{
    auto pc = PowerControlAmplifier.makeBlock(repeat(Complex!real(0, 0)), 10.dBm);
}