module dffdd.blockdiagram.amplifier;

import std.traits;
import std.range;
import std.math;
import std.complex;

import dffdd.utils.unit;
import dffdd.blockdiagram.utils;



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


struct PowerAmplifierConverterImpl(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;

    this(Gain gain, Voltage iip3, Voltage iip5 = Voltage(0)) pure nothrow @safe @nogc
    {
        _g1V = gain.gain;
        _g3V = iip3.V == 0 ? 0 : (gain.gain / iip3.V^^2);
        _g5V = iip5.V == 0 ? 0 : (gain.gain / iip5.V^^4);
    }


    void opCall(in C input, ref C output) const pure nothrow @safe @nogc
    {
        auto x = input;
        auto x1p = x.re^^2 + x.im^^2,
                x3 = x * x1p,
                x5 = x3 * x1p;

        output = x * _g1V + x3 * _g3V + x5 * _g5V;
    }


    // void opCall(R)(R input, ref C[] output) const pure nothrow @safe
    // if(isInputRange!R && hasLength!R && is(ElementType!R : C))
    // {
    //     immutable len = input.length;
    //     output.length = len;

    //     foreach(i; 0 .. len){
    //         this.opCall(input.front, output[i]);
    //         input.popFront();
    //     }
    // }


    PowerAmplifierConverterImpl dup() const pure nothrow @safe @nogc @property
    {
        return this;
    }


    private:
    real _g1V;
    real _g3V;
    real _g5V;
}


alias PowerAmplifierConverter(C) = ArrayConverterOf!(PowerAmplifierConverterImpl!C);


PowerAmplifierConverter!C makePowerAmplifier(C)(Gain gain, Voltage iip3, Voltage iip5 = Voltage(0))
{
    return PowerAmplifierConverter!C(gain, iip3, iip5);
}

unittest 
{
    import std.complex;
    import std.random;
    import std.algorithm;
    import dffdd.blockdiagram.utils;
    import std.math;
    import std.stdio : writeln;

    alias C = Complex!float;

    auto signal = new C[64];
    foreach(ref e; signal)
        e = C(uniform01(), uniform01);

    auto r0 = PowerAmplifier!(C[])(signal, 30.dB, 30.dBm, 30.dBm).array();
    auto r1 = signal.connectTo(makePowerAmplifier!C(30.dB, 30.dBm, 30.dBm));
    auto r2 = signal.chunks(3).connectTo(makePowerAmplifier!C(30.dB, 30.dBm, 30.dBm)).joiner;

    assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r1));
    assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r2));
}


/**
o: saturation value
s: smooth factor
u: |x[n]|
g(u) = (G*u)/(1+(u/o)^^(2s))^^(1/(2*s))
*/
struct RappModel(R)
{
    this(R r, Gain gain, real smoothFactor, real saturation)
    {
        _r = r;

        _g = gain.gain;
        _s = smoothFactor;
        _o = saturation;
    }


    auto front()
    {
        auto x = _r.front;
        auto r = abs(x),
             u = x / r;     // unit vector

        // rが小さすぎるときに，単位ベクトルが発散するのを防ぐ
        if(r <= 1E-6)
            return x;

        r = (r * _g) / (( 1 + (r/_o)^^(2*_s) )^^(1/(2*_s)));
        return r * u;
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
        dst._g = this._g;
        dst._s = this._s;
        dst._o = this._o;

        return dst;
    }
  }


  private:
    R _r;
    real _g, _s, _o;
}


struct RappModelConverterImpl(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(Gain gain, real smoothFactor, real saturation)
    {
        _g = gain.gain;
        _s = smoothFactor;
        _o = saturation;
    }


    void opCall(InputElementType input, ref OutputElementType output) const// pure nothrow @safe @nogc
    {
        auto x = input;
        auto r = sqAbs(x / _o);

        // rが小さすぎるときに，単位ベクトルが発散するのを防ぐ
        if(_s == 1){
            r = (_g) / (sqrt( 1 + r));
            output = r * x;
        }
        else if(_s == 3){
            r = (_g) / (sqrt( 1 + r^^3).cbrt);
            output = r * x;
        }
        else{
            r = (_g) / (( 1 + r^^(_s) )^^(1/(2*_s)));
            output = r * x;
        }
    }


    void opCall(R)(R input, ref OutputElementType[] output) const pure nothrow @trusted
    if(isInputRange!R && hasLength!R && is(ElementType!R : InputElementType))
    {
        output.length = input.length;

        foreach(i; 0 .. input.length){
            this.opCall(input.front, output[i]);
            input.popFront();
        }
    }


    RappModelConverterImpl dup() const pure nothrow @safe @nogc @property
    {
        return this;
    }


  private:
    real _g, _s, _o;
}


alias RappModelConverter(C) = ArrayConverterOf!(RappModelConverterImpl!C);


RappModelConverter!C makeRappModel(C)(Gain gain, real smoothFactor, real saturation)
{
    return RappModelConverter!C(gain, smoothFactor, saturation);
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

    foreach(s; 1 .. 5){
        auto r0 = RappModel!(C[])(signal, 30.dB, s, 1).array;
        auto r1 = signal.connectTo(makeRappModel!C(30.dB, s, 1));
        auto r2 = signal.chunks(3).connectTo(makeRappModel!C(30.dB, s, 1)).joiner;

        assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r1));
        assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r2));
    }
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


struct LinearAmplifierConverterImpl(C)
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


    LinearAmplifierConverterImpl dup() const @property pure nothrow @safe @nogc
    {
        return this;
    }


  private:
    real _gain1V;
}


alias LinearAmplifierConverter(C) = ArrayConverterOf!(LinearAmplifierConverterImpl!C);


LinearAmplifierConverter!C makeLinearAmplifier(C)(Gain gain)
{
    return LinearAmplifierConverter!C(gain);
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

    auto r0 = VGA!(C[])(signal, 30.dB).array;
    auto r1 = signal.connectTo(makeLinearAmplifier!C(30.dB));
    auto r2 = signal.chunks(2).connectTo(makeLinearAmplifier!C(30.dB)).joiner;

    assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r1));
    assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r2));
}


struct PowerControlAmplifier
{
    import std.stdio;

    static
    auto makeBlock(R)(R r, Voltage op, size_t avgSize = 512)
    {
        return PowerControlAmplifierImpl!R(r, op, avgSize);
    }


    static
    struct PowerControlAmplifierImpl(R)
    {
        this(R r, Voltage op, size_t avgSize)
        {
            _r = r;
            if(_r.empty){
                _empty = true;
                return;
            }

            _power = op.V^^2;
            _alpha = 1;
            _front = _r.front;
            _r.popFront();

            _empty = false;
            _cnt = 0;
            _avgSize = avgSize;
            _sumPower = 0;
            _avgCount = 0;
        }


        auto front() const @property
        {
            return _alpha * _front;
        }


        bool empty() const @property { return _empty; }


        void popFront()
        {
            if(_r.empty){
                _empty = true;
                return;
            }

            _sumPower += _front.re^^2 + _front.im^^2;
            ++_cnt;
            if(_cnt == _avgSize){
                if(_avgCount < 30){
                    if(_sumPower == 0)
                        _alpha = _alpha;
                    else
                        _alpha = _alpha / 2 + sqrt(_power / (_sumPower / _avgSize)) / 2;

                    ++_avgCount;
                }

                _sumPower = 0;
                _cnt = 0;
            }

            _front = this._r.front;

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
        real _power;
        real _alpha;
        Unqual!(ElementType!R) _front;
        bool _empty;
        size_t _cnt;
        size_t _avgSize;
        size_t _avgCount;
        real _sumPower;
    }
}

unittest
{
    auto pc = PowerControlAmplifier.makeBlock(repeat(Complex!real(0, 0)), 10.dBm);
}


struct PowerControlAmplifierConverterImpl(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(Voltage op, size_t avgSize)
    {
        _power = op.V^^2;
        _alpha = 1;

        _cnt = 0;
        _avgSize = avgSize;
        _sumPower = 0;
        _avgCount = 0;
    }


    void opCall(InputElementType input, ref OutputElementType output)// pure nothrow @safe @nogc
    {
        _sumPower += input.re^^2 + input.im^^2;
        ++_cnt;
        if(_cnt == _avgSize){
            if(_avgCount < 30){
                if(_sumPower == 0)
                    _alpha = _alpha;
                else
                    _alpha = _alpha / 2 + sqrt(_power / (_sumPower / _avgSize)) / 2;

                // writefln("%s : %s : %s", _power, _alpha, _sumPower);
                ++_avgCount;
            }

            _sumPower = 0;
            _cnt = 0;
        }

        output = input * _alpha;
    }


    PowerControlAmplifierConverterImpl dup() const // pure nothrow @safe @nogc @property
    {
        return this;
    }


  private:
    real _power;
    real _alpha;
    size_t _cnt;
    size_t _avgSize;
    size_t _avgCount;
    real _sumPower;
}


alias PowerControlAmplifierConverter(C) = ArrayConverterOf!(PowerControlAmplifierConverterImpl!C);


PowerControlAmplifierConverter!C makePowerControlAmplifier(C)(Voltage op, size_t avgSize = 512)
{
    return PowerControlAmplifierConverter!C(op, avgSize);
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

    auto r0 = PowerControlAmplifier.makeBlock(signal, 30.dBm, 512).array;
    auto r1 = signal.connectTo(makePowerControlAmplifier!C(30.dBm, 512));
    auto r2 = signal.chunks(3).connectTo(makePowerControlAmplifier!C(30.dBm, 512)).joiner;

    assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r1));
    assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r2));
}
