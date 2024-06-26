module dffdd.blockdiagram.amplifier;

import std.algorithm;
import std.traits;
import std.range;
import std.math;
import std.complex;
import std.json;

import dffdd.blockdiagram.utils;
import dffdd.utils.unit;
import dffdd.utils.json;
import dffdd.math.math;


interface IAmplifier(C) : IConverter!(C, C, Yes.isDuplicatable)
{
    Gain linearGain() const;
    JSONValue dumpInfoToJSON() const;
}


final class AmplifierObject(C, Struct) : IAmplifier!C
{
    this(Struct original)
    {
        _original = original;
    }


    alias InputElementType = Struct.InputElementType;
    alias OutputElementType = Struct.OutputElementType;


    void opCall(in InputElementType[] src, ref OutputElementType[] dst)
    {
        static if(is(typeof(_original(src, dst))))
            _original(src, dst);
        else {
            foreach(i; 0 .. src.length) {
                _original(src[i], dst[i]);
            }
        }
    }


    AmplifierObject!(C, Struct) dup() const
    {
        return new AmplifierObject!(C, Struct)(_original.dup);
    }


    Gain linearGain() const
    {
        return _original.linearGain();
    }


    JSONValue dumpInfoToJSON() const
    {
        return _original.dumpInfoToJSON();
    }


  private:
    Struct _original;
}


IAmplifier!(C) toAmplifierObject(C, Amp)(Amp amp)
{
    return new AmplifierObject!(C, Amp)(amp);
}


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
Rapp model
*/
struct RappModelConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;
    alias F = typeof(C.init.re);

    /**
    gain: 小信号ゲイン
    smoothFactor: スムースネスファクタ
    outputSatV: 出力飽和電圧（電力）
    */
    this(F smoothFactor, Gain gain, Voltage outputSatV)
    {
        _g = gain.asV;
        _s = smoothFactor;
        _o = outputSatV.volt;

        if(_s == cast(int)_s) {
            _bInteger = true;
            _int_s = cast(int)_s;
        }
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        auto r = fast_abs(input),
             u = input / r;         // unit vector

        // rが小さすぎるときに，単位ベクトルが発散するのを防ぐ
        if(r <= 1E-6){
            output = input * _g;
        }else{
            r *= _g;

            F p;    // (r/_o)^^(2*s)
            if(_bInteger) {
                p = fast_powi!F(r / _o, 2 * _int_s);
            } else {
                p = fast_pow!F(r / _o, 2 * _s);
            }

            F q;    // (1 + p)^^(1/(2*s))
            if(_bInteger && _int_s == 1) {
                q = fast_sqrt!F(1 + p);
            } else if(_bInteger && _int_s == 2) {
                q = fast_sqrt!F(fast_sqrt!F(1 + p));
            } else if(_bInteger && _int_s == 3) {
                q = fast_cbrt!F(fast_sqrt!F(1 + p));
            } else if(_bInteger && _int_s == 4) {
                q = fast_sqrt!F(fast_sqrt!F(fast_sqrt!F(1 + p)));
            } else {
                q = fast_pow!F(1 + p, 1.0L / (2 * _s));
            }

            r /= q;
            output = r * u;
        }
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        if(output.length != input.length)
            output.length = input.length;

        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    /**
    線形領域での利得を返します
    */
    Gain linearGain() const @property
    {
        return Gain.fromVoltageGain(_g);
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
            "outputSaturationVoltage":   _o
        ]);
    }


  private:
    F _g, _s, _o;
    bool _bInteger;
    int _int_s;
}


struct RappModel
{
    static
    auto makeBlock(R)(R range, real smoothFactor, Gain gain, Voltage outputSatV)
    if(isInputRange!R)
    {
        import dffdd.blockdiagram.utils : connectTo;
        alias E = Unqual!(ElementType!R);
        return range.connectTo!(RappModelConverter!E)(smoothFactor, gain, outputSatV);
    }


    static
    auto makeBlock(C)(BufferedOutputTerminal!C src, size_t len, real smoothFactor, Gain gain, Voltage outputSatV)
    {
        auto conv = RappModelConverter!C(smoothFactor, gain, outputSatV);
        return new ConverterBlock!RappModelConverter(len, src, conv);
    }
}

unittest
{
    auto r = RappModel.makeBlock(Complex!float[].init, 1, 1.0.dB, 30.dBm);
}



/**
Salehモデルを構築します．
Salehモデルは，4つのパラメータAa, Ba, Ap, Bpで以下のように表されます．

f(x) = Aa x / (1 + Ba |x|^2) * expi ( Ap |x|^2 / (1 + Bp |x|^2) )

この関数の振幅は，|x|^2 = 1/Ba で最大値である f(1 / sqrt{Ba}) = Aa/(2 sqrt(Ba)) になります．
また，小信号ゲインは Aa です．

この実装では，この振幅の最大値と小信号ゲインをユーザーが変更可能にしています．
*/
struct SalehModelConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;

    /**
    gain: 小信号ゲイン
    osatV: 出力飽和電圧
    satPhi: 出力飽和電圧での位相回転量
    */
    this(Gain gain, Voltage osatV, real satPhi)
    {
        _g = gain.asV;
        _o = osatV.volt;
        _phi = satPhi;
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        immutable r = std.complex.abs(input),
                  u = input / r;    // unit vector

        // rが小さすぎるときに，単位ベクトルが発散するのを防ぐ
        if(r <= 1E-6) {
            output = input * _g;
        } else {
            output = _o * normalized_saleh(_g * r / _o, _phi) * u;
        }
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        if(output.length != input.length)
            output.length = input.length;

        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    /**
    線形領域での利得を返します
    */
    Gain linearGain() const @property
    {
        return Gain.fromVoltageGain(_g);
    }


    typeof(this) dup() const pure nothrow @safe @nogc @property
    {
        return this;
    }


    JSONValue dumpInfoToJSON() const
    {
        return JSONValue([
            "gain":         _g,
            "saturation":   _o,
            "phisat" : _phi
        ]);
    }


  private:
    real _g, _o, _phi;


    /**
    飽和電圧1，小信号ゲイン1，飽和時の位相回転量phiのsalehモデル
    */
    C normalized_saleh(real r, real phi)
    {
        immutable real aa = 1.0,
                       ba = aa^^2 / 4,
                       ap = 0.5 * phi,
                       bp = ba;

        return C(aa * r / (1 + ba*r^^2) * std.complex.expi(ap*r^^2 / (1+bp*r^^2)));
    }
}


/**
o: input saturation value
*/
struct SoftLimitConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;

    this(Gain gain, Voltage osatV)
    {
        _g = gain.asV;
        _o = osatV.volt;
    }


    void opCall(InputElementType input, ref OutputElementType output) @nogc
    {
        immutable r = std.complex.abs(input),
                  u = input / r;    // unit vector

        if(_g * r < _o) {
            output = _g * input;
        } else {
            output = _o * u;
        }
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        if(output.length != input.length)
            output.length = input.length;

        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    /**
    線形領域での利得を返します
    */
    Gain linearGain() const @property
    {
        return Gain.fromVoltageGain(_g);
    }


    typeof(this) dup() const pure nothrow @safe @nogc @property
    {
        return this;
    }


    JSONValue dumpInfoToJSON() const
    {
        return JSONValue([
            "gain":         _g,
            "saturation":   _o
        ]);
    }


  private:
    real _g, _o;
}



struct VGAConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;

    this(Gain gain)
    {
        _gain1V = gain.asV;
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


    /**
    64サンプルからスタートして,
    64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 32768, 32768, ...サンプルごとに平均電力値の計算をし，増幅率を更新します
    前回，前々回からのゲインの増減率が1±rtol以下であれば停止します
    */
    this(Voltage op, real rtol = 1E-2, size_t startSamples = 64, size_t totalCount = 10)
    {
        _targetPower = op.volt^^2;
        _gain = 1;
        _lastGain = 0;
        _cnt = 0;
        _rtol = rtol;
        _sumPower = 0;
        _avgTotal = startSamples;
        _maxTotal = 64 << totalCount;
    }


    void opCall(InputElementType input, ref OutputElementType output) @nogc
    {
        if(!_isConverged){
            _sumPower += input.re^^2 + input.im^^2;
            ++_cnt;
            if(_cnt == _avgTotal) {
                if(_sumPower != 0){
                    real lastGain = _lastGain;
                    real nowGain = _gain;
                    real newGain = sqrt(_targetPower / (_sumPower / _avgTotal));

                    if(abs(lastGain / newGain - 1) < _rtol
                    && abs(nowGain / newGain - 1) < _rtol)
                    {
                        _isConverged = true;
                    }

                    _gain = newGain;
                    _lastGain = nowGain;
                }
                
                if(_avgTotal < _maxTotal)
                    _avgTotal *= 2;

                _cnt = 0;
                _sumPower = 0;
            }
        }

        output = input * _gain;
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        if(output.length != input.length)
            output.length = input.length;

        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    Gain gain() @property const
    {
        return Gain.fromVoltageGain(_gain);
    }


    bool isConverged() @property const
    {
        return _isConverged;
    }


    PowerControlAmplifierConverter!C dup() const pure nothrow @safe @nogc @property 
    {
        return this;
    }


    JSONValue dumpInfoToJSON() const
    {
        return JSONValue([
            "gain": _gain
        ]);
    }


  private:
    real _targetPower;
    real _gain;
    real _rtol;
    real _lastGain;
    real _sumPower;
    size_t _cnt;
    size_t _avgTotal;
    size_t _maxTotal;
    bool _isConverged;
}



struct PowerControlAmplifier
{
    static
    auto makeBlock(R)(R r, Voltage op, real rtol = 1E-2, size_t startSamples = 64, size_t totalCount = 10)
    {
        import dffdd.blockdiagram.utils : connectTo;
        alias E = Unqual!(ElementType!R);
        return r.connectTo!(PowerControlAmplifierConverter!E)(op, rtol, startSamples, totalCount);
    }
}

unittest
{
    auto pc = PowerControlAmplifier.makeBlock(repeat(Complex!real(0, 0)), 10.dBm);
}



/**
ret[0]: input gain
ret[1]: output gain
*/
private R[2] convertLinGainAndVsat(R)(R linGain, R Vsat, R targetLinGain, R targetVsat)
{
    return cast(R[2])[
        Vsat / targetVsat * targetLinGain / linGain,
        targetVsat / Vsat
    ];
}

unittest
{
    alias C = Complex!double;

    auto sf = SoftLimitConverter!C(0.dB, 30.dBm);
    auto gs = convertLinGainAndVsat(0.dB.asV, 30.dBm.volt, 20.dB.asV, 0.dBm.volt);

    C o;
    sf(C(100), o);
    assert(o == C(1));
    sf(C(0.01), o);
    assert(o == C(0.01));

    C newx = C(100) * gs[0];
    sf(newx, o);
    o *= gs[1];
    assert(o.re.isClose(C(0.001.sqrt)));
    assert(o.im.isClose(0));

    newx = C(0.0001) * gs[0];
    sf(newx, o);
    o *= gs[1];
    assert(o.re.isClose(0.001));
    assert(o.im.isClose(0));
}



/**
与えられたAM/AM特性を線形補間する増幅器モデル
*/
struct LinearInterpolatedAMAMConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;
    alias F = typeof(C.init.re);


    /**
    入力信号の平均電力が1のときに最適であるとして設計されたAM/AM特性カーブ（xs, fs）を線形補間する増幅器を作成します．
    また，このカーブに基づいて，小信号ゲインと飽和電圧をユーザーが変更可能にしています．
    */
    this(immutable(F)[] xs, immutable(F)[] fs, Gain gain, Voltage osatV)
    {
        _g = gain.asV;
        _o = osatV.volt;
        _xs = xs;
        _fs = fs;
        _interp = _makeInterpolationObject(xs, fs);
        _x0 = xs[0];
        _f0 = fs[0];
        _xN = xs[$-1];
        _fN = fs[$-1];
        _maxf = fs.maxElement;

        // 正規化前の小信号ゲインの計算
        if(0.01 <= _x0) {
            _g0 = _f0 / _x0;
        } else {
            _g0 = _interp(0.01) / 0.01;
        }
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        F inputAmp = std.complex.abs(input);
        C u = input / inputAmp;
        if(u.re.isNaN || u.im.isNaN) {
            output = input * _g;
        } else {
            output = normalized_AMAM(inputAmp * _g / _o) * _o * u;
        }
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        if(output.length != input.length)
            output.length = input.length;

        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    /**
    線形領域での利得を返します
    */
    Gain linearGain() const @property
    {
        return Gain.fromVoltageGain(_g);
    }


    LinearInterpolatedAMAMConverter!C dup() const
    {
        return LinearInterpolatedAMAMConverter!C(_xs, _fs, Gain.fromVoltageGain(_g), Voltage(_o));
    }


    JSONValue dumpInfoToJSON() const
    {
        JSONValue jv = JSONValue(["gain": _g, "outputSaturationVoltage": _o]);
        jv["xs"] = _xs;
        jv["fs"] = _fs;
        return jv;
    }


  private:
    alias InterpT = ReturnType!_makeInterpolationObject;
    InterpT _interp;
    F _g, _o;
    immutable(F)[] _xs, _fs;
    F _x0, _f0, _xN, _fN, _g0;
    F _maxf;

    static
    auto _makeInterpolationObject(immutable(F)[] xs, immutable(F)[] fs)
    {
        import mir.interpolate.linear : linear;
        import mir.rc.array;
        return linear!F(xs.rcarray.asSlice, fs.rcarray!(const(F)).asSlice);
    }


    /**
    小信号ゲインを 1 ，出力飽和電圧を 1 にしたAM/AM特性
    */
    real normalized_AMAM(real r)
    {
        immutable r_normalized = r * _maxf / _g0;

        if(r_normalized <= _x0) {
            return r;
        } else if(r_normalized >= _xN) {
            return _fN / _maxf;
        } else {
            return _interp(r_normalized) / _maxf;
        }
    }
}

unittest
{
    alias F = double;
    alias C = Complex!double;

    immutable(F)[] xs = [0.10, 0.5, 1.0];
    immutable(F)[] fs = [0.04, 0.2, 0.8];

    Complex!F y;

    // f(x) = 0.25 x        when x = [0, 4)
    // f(x) = 0.75 x - 2    when x = [4, 8)
    // f(x) = 4             when x = [8, inf)
    auto interp = LinearInterpolatedAMAMConverter!C(xs, fs, Gain.fromVoltageGain(0.25), Voltage(4));

    interp(C(0), y);
    assert(y.re.isClose(0));

    interp(C(2), y);
    assert(y.re.isClose(0.5));

    interp(C(4), y);
    assert(y.re.isClose(1));
    
    interp(C(6), y);
    assert(y.re.isClose(3*1.5 - 2));

    interp(C(8), y);
    assert(y.re.isClose(4));

    interp(C(9), y);
    assert(y.re.isClose(4));

    interp(C(10), y);
    assert(y.re.isClose(4));
}


/**
多項式で近似された増幅器モデル
*/
struct PolynomialWithClipConverter(C)
{
    private {
        alias R = typeof(C.init.re);
    }

    alias InputElementType = C;
    alias OutputElementType = C;


    this(immutable(C)[] cs, Gain gain, Voltage osatV)
    in(cs.length >= 1)
    {
        _cs = cs;
        _targetGain = gain.asV;
        _clipLevel = osatV.volt;

        _convGains = convertLinGainAndVsat!R(_cs[0].abs, 1, gain.asV, osatV.volt);
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        R inputAmp = std.complex.abs(input);
        C u = input / inputAmp;
        if(u.re.isNaN || u.im.isNaN) {
            output = input * _targetGain;
        } else {
            import dffdd.math : poly;

            immutable r = inputAmp * _convGains[0];
            output = poly(C(r^^2), _cs) * r * _convGains[1] * u;
        }

        if(output.sqAbs > _clipLevel^^2) {
            output = _clipLevel * u;
        }
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        if(output.length != input.length)
            output.length = input.length;

        foreach(i, e; input)
            this.opCall(e, output[i]);
    }


    Gain linearGain() const
    {
        return Gain.fromVoltageGain(_targetGain);
    }


    JSONValue dumpInfoToJSON() const
    {
        import std.base64 : Base64;
        import msgpack;

        return JSONValue([
            "cs_msgpackbin_base64": JSONValue(Base64.encode(msgpack.pack(_cs))),
            "gain": JSONValue(_targetGain),
            "clipLevel": JSONValue(_clipLevel),
        ]);
    }


    PolynomialWithClipConverter!C dup() const
    {
        return PolynomialWithClipConverter!C(_cs, Gain.fromVoltageGain(_targetGain), Voltage(_clipLevel));
    }


  private:
    immutable(C)[] _cs;

    R _targetGain;
    R _clipLevel;
    R[2] _convGains;
}

unittest
{
    alias R = double;
    alias C = Complex!R;
    

    // linear
    immutable(C)[] cs = [C(2)];

    auto conv = PolynomialWithClipConverter!C(cs, 0.dB, 30.dBm);

    C inp, outp;
    inp = C(0.1);
    conv(inp, outp);
    assert(outp.re.isClose(0.1));
    assert(outp.im.isClose(0));

    inp = C(100);
    conv(inp, outp);
    assert(outp.re.isClose(1));
    assert(outp.im.isClose(0));
}

