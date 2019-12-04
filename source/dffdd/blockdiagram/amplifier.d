module dffdd.blockdiagram.amplifier;

import std.algorithm;
import std.traits;
import std.range;
import std.math;
import std.complex;
import std.json;

import dffdd.utils.unit;
import dffdd.utils.json;
import dffdd.math.math;

import dffdd.blockdiagram.utils;


interface IAmplifier(C) : ISISOConverter!(C, C, Yes.isDuplicatable)
{
    Gain linearGain() const;
    JSONValue dumpInfoToJSON() const;
}


/**
Rapp model
*/
struct RappModel(C)
{
    alias IType = C;
    alias OType = C;
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


    void opCall(IType input, ref OType output)
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


    void opCall(in IType[] input, ref OType[] output)
    {
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


/**
Salehモデルを構築します．
Salehモデルは，4つのパラメータAa, Ba, Ap, Bpで以下のように表されます．

f(x) = Aa x / (1 + Ba |x|^2) * expi ( Ap |x|^2 / (1 + Bp |x|^2) )

この関数の振幅は，|x|^2 = 1/Ba で最大値である f(1 / sqrt{Ba}) = Aa/(2 sqrt(Ba)) になります．
また，小信号ゲインは Aa です．

この実装では，この振幅の最大値と小信号ゲインをユーザーが変更可能にしています．
*/
struct SalehModel(C)
{
    alias IType = C;
    alias OType = C;

    /**
    gain: 小信号ゲイン
    osatV: 出力飽和電圧
    */
    this(Gain gain, Voltage osatV)
    {
        _g = gain.asV;
        _o = osatV.volt;
    }


    void opCall(IType input, ref OType output)
    {
        input *= _g;

        immutable r = std.complex.abs(input),
                  u = input / r;    // unit vector

        // rが小さすぎるときに，単位ベクトルが発散するのを防ぐ
        if(r <= 1E-6) {
            output = input * _g;
        } else {
            output = _o * normalized_saleh(_g * r / _o) * u;
        }
    }


    void opCall(in IType[] input, ref OType[] output)
    {
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


    /**
    飽和電圧1，小信号ゲイン1のsalehモデル
    */
    C normalized_saleh(real r)
    {
        immutable aa = 2.1587,
                  ba = aa^^2 / 4,   // = 1.16402521
                  ap = 4.033,
                  bp = 9.1040;

        r /= aa;
        return C(aa * r / (1+ba*r^^2) * std.complex.expi(ap*r^^2/(1+bp*r^^2)));
    }
}


/**
o: input saturation value
*/
struct IdealAmplifier(C)
{
    alias IType = C;
    alias OType = C;

    this(Gain gain, Voltage osatV)
    {
        _g = gain.asV;
        _o = osatV.volt;
    }


    void opCall(IType input, ref OType output)
    {
        immutable r = std.complex.abs(input),
                  u = input / r;    // unit vector

        if(_g * r < _o) {
            output = _g * input;
        } else {
            output = _o * u;
        }
    }


    void opCall(in IType[] input, ref OType[] output)
    {
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



struct FixedGainAmplifier(C)
{
    alias IType = C;
    alias OType = C;

    this(Gain gain)
    {
        _gain1V = gain.asV;
    }


    void opCall(IType input, ref OType output)
    {
        output = input * _gain1V;
    }


    void opCall(in IType[] input, ref OType[] output)
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


    FixedGainAmplifier!C dup() const pure nothrow @safe @nogc @property 
    {
        return this;
    }


  private:
    real _gain1V;
}

unittest
{
    alias T = FixedGainAmplifier!float;
    static assert(isConverter!T);
    static assert(isSISOConverter!T);
}


struct PowerControlAmplifier(C)
{
    alias IType = C;
    alias OType = C;


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


    void opCall(IType input, ref OType output) @nogc
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


    void opCall(in IType[] input, ref OType[] output)
    {
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


    PowerControlAmplifier!C dup() const pure nothrow @safe @nogc @property 
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


/**
与えられたAM/AM特性を線形補間する増幅器モデル
*/
struct LinearInterpolatedAMAM(C)
{
    alias IType = C;
    alias OType = C;
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


    void opCall(IType input, ref OType output)
    {
        F inputAmp = std.complex.abs(input);
        C u = input / inputAmp;
        if(u.re.isNaN || u.im.isNaN) {
            output = input * _g;
        } else {
            output = normalized_AMAM(inputAmp * _g / _o) * _o * u;
        }
    }


    void opCall(in IType[] input, ref OType[] output)
    {
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


    LinearInterpolatedAMAM!C dup() const
    {
        return LinearInterpolatedAMAM!C(_xs, _fs, Gain.fromVoltageGain(_g), Voltage(_o));
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
        import mir.ndslice;
        return linear!F(xs.sliced, fs.sliced);
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
    auto interp = LinearInterpolatedAMAM!C(xs, fs, Gain.fromVoltageGain(0.25), Voltage(4));

    interp(C(0), y);
    assert(y.re.approxEqual(0));

    interp(C(2), y);
    assert(y.re.approxEqual(0.5));

    interp(C(4), y);
    assert(y.re.approxEqual(1));
    
    interp(C(6), y);
    assert(y.re.approxEqual(3*1.5 - 2));

    interp(C(8), y);
    assert(y.re.approxEqual(4));

    interp(C(9), y);
    assert(y.re.approxEqual(4));

    interp(C(10), y);
    assert(y.re.approxEqual(4));
}
