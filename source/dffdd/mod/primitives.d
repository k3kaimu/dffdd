module dffdd.mod.primitives;

import std.complex;
import std.math;
import std.numeric : gcd;
import std.traits;
import std.typecons;

import dffdd.math.math;



enum bool isModulator(T) = is(typeof((T t){
    alias IType = T.InputElementType;
    alias OType = T.OutputElementType;

    size_t ips = t.symInputLength;      // 1シンボル出力するために必要な入力のサイズ
    size_t ops = t.symOutputLength;     // 1シンボル出力したときの出力のサイズ

    IType[] inputs;
    OType[] outputs;

    outputs = t.modulate(cast(const)inputs, outputs);
    inputs = t.demodulate(cast(const)outputs, inputs);
}));


Mod.OutputElementType[] mod(Mod, E)(ref Mod mod, in E[] sig)
{
    typeof(return) dst;

    mod.modulate(sig, dst);

    return dst;
}


Mod.InputElementType[] demod(Mod, E)(ref Mod mod, in E[] sig)
{
    typeof(return) dst;

    mod.demodulate(sig, dst);

    return dst;
}


struct Bit
{
    ubyte value;
    alias value this;


    this(long v)
    {
        this.opAssign(v);
    }


    void opAssign(long v)
    {
        if(v)
            value = 1;
        else
            value = 0;
    }
}


struct ChainedMod(Mod1, Mod2)
if(isModulator!Mod1 && isModulator!Mod2 && is(Mod1.OutputElementType == Mod2.InputElementType))
{
    alias InputElementType = Mod1.InputElementType;
    alias OutputElementType = Mod2.OutputElementType;


    this(Mod1 mod1, Mod2 mod2)
    {
        _mod1 = mod1;
        _mod2 = mod2;
    }


    size_t symInputLength() const @property { return _mod1.symInputLength * lcm(_mod1.symOutputLength, _mod2.symInputLength) / _mod1.symOutputLength; }
    size_t symOutputLength() const @property { return _mod2.symOutputLength * lcm(_mod1.symOutputLength, _mod2.symInputLength) / _mod2.symInputLength; }


    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    in{
        assert(inputs.length % this.symInputLength == 0);
        assert(outputs.length % this.symOutputLength == 0);
    }
    do{
        Mod1.OutputElementType[] buf;
        _mod1.modulate(inputs, buf);
        _mod2.modulate(buf, outputs);

        return outputs;
    }


    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    in{
        assert(outputs.length % this.symInputLength == 0);
        assert(inputs.length % this.symOutputLength == 0);
    }
    do{
        Mod2.InputElementType[] buf;
        _mod2.demodulate(inputs, buf);
        _mod1.demodulate(buf, outputs);

        return outputs;
    }


  private:
    Mod1 _mod1;
    Mod2 _mod2;


    static
    T lcm(T)(T a, T b)
    {
        import std.numeric : gcd;
        return a * (b / gcd(a, b));
    }
}


auto chainedMod(Mod1, Mod2)(Mod1 mod1, Mod2 mod2)
{
    return ChainedMod!(Mod1, Mod2)(mod1, mod2);
}


///
unittest
{
    import std.algorithm : map;
    import std.array : array;
    import std.complex;
    import dffdd.mod.ofdm;
    import dffdd.mod.qpsk;

    // 1次変調: QPSK
    // FFTサイズ: 8
    // サイクリックプレフィックス: 3
    // 使用サブキャリア数: 4
    // アップサンプリング率: 2
    auto mod = chainedMod(new QPSK!(Complex!float)(), new OFDM!(Complex!float)(8, 3, 4, 2));

    assert(mod.symInputLength == 2 * 4);
    assert(mod.symOutputLength == (8 + 3) * 2);

    Bit[] inps = [1, 1, 0, 0, 1, 0, 0, 1].map!(a => Bit(a & 0xFF)).array();
    Complex!float[] res;

    // 変調
    mod.modulate(inps, res);

    // CPのチェック
    foreach(i; 0 .. 3 * 2)
        assert(res[i] == res[$ - 3*2 + i]);

    Bit[] decoded;

    // 復調
    mod.demodulate(res, decoded);

    // 復調結果のチェック
    assert(decoded == inps);
}



/**
変復調機を変調器か復調器のConverterにします
*/
struct ModConverter(Mod, Flag!"isMod" isMod)
{
    this(Mod mod)
    {
        _mod = mod;
    }


    static if(isMod)
    {
        alias InputElementType = Mod.InputElementType;
        alias OutputElementType = Mod.OutputElementType;
    }
    else
    {
        alias InputElementType = Mod.OutputElementType;
        alias OutputElementType = Mod.InputElementType;
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        static if(isMod)
        {
            _mod.modulate(input, output);
        }
        else
        {
            _mod.demodulate(input, output);
        }
    }


  private:
    Mod _mod;
}


/// ditto
ModConverter!(Mod, Yes.isMod) asModConverter(Mod)(Mod mod)
{
    return ModConverter!(Mod, Yes.isMod)(mod);
}


/// ditto
ModConverter!(Mod, No.isMod) asDemodConverter(Mod)(Mod mod)
{
    return ModConverter!(Mod, No.isMod)(mod);
}


unittest 
{
    static struct DummyMod(R)
    {
        alias InputElementType = R;
        alias OutputElementType = R;

        void modulate(in R[] input, ref R[] output) {
            output.length = input.length;
            output[] = input[] * 2;
        }


        void demodulate(in R[] input , ref R[] output) {
            output.length = input.length;
            output[] = input[] / 2;
        }
    }

    auto dummy = DummyMod!float();
    auto mod = dummy.asModConverter;
    auto demod = dummy.asDemodConverter;

    float[] input = [1, 2, 3];
    float[] output;
    mod(input, output);
    assert(output == [2, 4, 6]);

    demod(input, output);
    assert(output == [0.5, 1, 1.5]);
}


size_t countBitError(in ushort[] symAs, in ushort[] symBs)
in(symAs.length == symBs.length)
{
    import core.bitop;

    size_t cnt = 0;
    foreach(i; 0 .. symAs.length)
        cnt += popcnt(symAs[i] ^ symBs[i]);

    return cnt;
}

///
unittest 
{
    ushort[] as = [0b0011, 0b1100, 0b0101, 0b1111];
    ushort[] bs = [0b0011, 0b1111, 0b1010, 0b1110];

    assert(countBitError(as, as) == 0);
    assert(countBitError(as, bs) == 7);
}



struct BERCounter
{
    this(size_t k)
    {
        _k = k;
    }


    void count(in ushort[] symAs, in ushort[] symBs)
    in(symAs.length == symBs.length)
    {
        _totalSym += symAs.length;
        _errBit += countBitError(symAs, symBs);
        foreach(i; 0 .. symAs.length)
            _errSym += (symAs[i] != symBs[i] ? 1 : 0);
    }


    struct BERCounterResult
    {
        size_t totalBits;
        size_t totalSyms;
        size_t errBits;
        size_t errSyms;
        double ber;
        double ser;
    }


    BERCounterResult result() const
    {
        BERCounterResult res;

        res.totalBits = _totalSym * _k;
        res.totalSyms = _totalSym;
        res.errBits = _errBit;
        res.errSyms = _errSym;
        res.ber = 1.0 * _errBit / (_totalSym * _k);
        res.ser = 1.0 * _errSym / _totalSym;

        return res;
    }


    size_t totalSyms() const { return _totalSym; }
    size_t totalBits() const { return _totalSym * _k; }
    size_t errBits() const { return _errBit; }
    size_t errSyms() const { return _errSym; }


  private:
    size_t _k;
    size_t _totalSym;
    size_t _errSym;
    size_t _errBit;
}

unittest
{
    ushort[] as = [0b0011, 0b1100, 0b0101, 0b1111];
    ushort[] bs = [0b0011, 0b1111, 0b1010, 0b1110];

    BERCounter counter = BERCounter(4); // one symbol consists of four bits.

    counter.count(as, bs);

    with(counter.result) {
        assert(totalBits == 16);
        assert(totalSyms == 4);
        assert(errBits == 7);
        assert(errSyms == 3);
        assert(ber == 7.0 / 16);
        assert(ser == 3.0 / 4);
    }
}



/**変調器の全Constellationを返します．mod.symOutputLength == 1でないといけません
*/
Mod.OutputElementType[] allConstellationPoints(Mod)(Mod mod)
if(isModulator!Mod)
in(mod.symOutputLength == 1)
{
    alias C = Mod.OutputElementType;

    size_t nbits = mod.symInputLength;

    C[] sym = new C[1];
    Bit[] bs = new Bit[nbits];
    C[] points = [];
    foreach(n; 0 .. 1 << nbits) {
        foreach(k; 0 .. nbits)
            bs[k] = (n & (1 << k)) ? 1 : 0;

        mod.modulate(bs, sym);
        points ~= sym[0];
    }

    return points;
}

unittest
{
    import dffdd.mod.bpsk;
    import dffdd.mod.qpsk;
    import dffdd.mod.qam;

    alias C = Complex!float;

    assert(BPSK!C().allConstellationPoints == [C(1, 0), C(-1, 0)]);
    assert(QPSK!C().allConstellationPoints
        == [C(SQRT1_2, SQRT1_2), C(-SQRT1_2, SQRT1_2), C(SQRT1_2, -SQRT1_2), C(-SQRT1_2, -SQRT1_2)]);
}



struct P0P1Calculator(C, F)
{
    this(Mod)(Mod mod)
    in(mod.symOutputLength == 1)
    {
        _nbits = mod.symInputLength;
        _p0 = new F[_nbits];
        _p1 = new F[_nbits];
        _points = mod.allConstellationPoints();
    }


    ref F[] computeP0P1(in C[] inputs, return ref F[] p0p1, F N0, C hcoef = C(1))
    {
        p0p1.length = inputs.length * _nbits;

        foreach(i, x; inputs) {
            _p0[] = 0;
            _p1[] = 0;

            foreach(j, y; _points) {
                // import dffdd.math.approxmath;
                F p = fast_exp!F(-(x - y * hcoef).sqAbs/ N0);

                foreach(k; 0 .. _nbits) {
                    if(j & (1 << k))
                        _p1[k] += p;
                    else
                        _p0[k] += p;
                }
            }

            foreach(k; 0 .. _nbits)
                p0p1[i * _nbits + k] = _p0[k] / _p1[k];
        }

        return p0p1;
    }


  private:
    size_t _nbits;
    C[] _points;
    F[] _p0;
    F[] _p1;
}


auto p0p1Calculator(F, Mod)(Mod mod)
{
    return P0P1Calculator!(Mod.OutputElementType, F)(mod);
}



struct SoftSymbolGenerator(C, R)
{
    this(Mod)(Mod mod)
    in(mod.symOutputLength == 1)
    do {
        _nbits = mod.symInputLength;
        _symbols = new C[2^^_nbits];
        _tempP0 = new R[_nbits];
        _tempP1 = new R[_nbits];
        
        Bit[] bits;
        foreach(n; 0 .. 2^^_nbits) {
            foreach(j; 0 .. _nbits) {
                bits ~= Bit(n & (1 << j));
            }
        }

        _symbols = mod.modulate(bits, _symbols);
    }


    private static R[] dontStorePower;


    ref C[] opCall(in R[] p0p1s, return ref C[] outputs, ref R[] powers = dontStorePower)
    in(p0p1s.length % _nbits == 0)
    do {
        outputs.length = p0p1s.length / _nbits;
        if(powers !is dontStorePower) {
            powers.length = p0p1s.length / _nbits;
            powers[] = 0;
        }

        foreach(i; 0 .. outputs.length) {
            outputs[i] = C(0);

            foreach(j; 0 .. _nbits) {
                _tempP1[j] = 1/(p0p1s[i * _nbits + j] + 1);
                _tempP0[j] = 1 - _tempP1[j];
            }

            foreach(n; 0 .. 2^^_nbits) {
                R q = 1;
                foreach(j; 0 .. _nbits) {
                    if(n & (1 << j))
                        q *= _tempP1[j];
                    else
                        q *= _tempP0[j];
                }

                outputs[i] += q * _symbols[n];

                if(powers !is dontStorePower)
                    powers[i] += q * _symbols[n].sqAbs;
            }
        }

        return outputs;
    }

  private:
    size_t _nbits;
    C[] _symbols;
    R[] _tempP0, _tempP1;
}


auto softSymbolGenerator(R, Mod)(Mod mod)
{
    return SoftSymbolGenerator!(Mod.OutputElementType, R)(mod);
}



struct LR01(F, alias mapping = [1, -1])
if(isFloatingPoint!F)
{
    F value;


    static
    LR01!(F, mapping) fromLLR(F llr)
    {
        return LR01!(F, mapping)(exp(llr));
    }


    // compute soft bit
    F softbit() const
    {
        auto p0 = (value == F.infinity) ? 1 : value / (value + 1);
        auto p1 = 1 - p0;
        return p0 * mapping[0] + p1 * mapping[1];
    }


    LLR01!(F, mapping) toLLR() const
    {
        return LLR01!(F, mapping)(log(value));
    }


    F llr() const
    {
        return this.toLLR().value;
    }


    typeof(this) toLR() const
    {
        return this;
    }


    F lr() const
    {
        return value;
    }


    int opCmp(G)(in LR01!(G, mapping) rhs) const
    {
        if(this.value == rhs.value)
            return 0;
        else if(this.value > rhs.value)
            return 1;
        else
            return -1;
    }


    LR01!(CommonType!(F, G), mapping) opBinary(string op : "*", G)(in LR01!(G, mapping) rhs) const
    {
        return typeof(return)(this.value * rhs.value);
    }


    LR01!(CommonType!(F, G), mapping) opBinary(string op : "/", G)(in LR01!(G, mapping) rhs) const
    {
        return typeof(return)(this.value / rhs.value);
    }


    LR01!(CommonType!(F, G), mapping) opBinary(string op : "&", G)(in LR01!(G, mapping) rhs) const
    {
        return this.opBinary!"*"(rhs);
    }


    LR01!(CommonType!(F, G), mapping) opBinary(string op : "|", G)(in LR01!(G, mapping) rhs) const
    {
        return this.opBinary!"/"(rhs);
    }
}

unittest
{
    auto x1 = LR01!float(float.infinity);
    assert(x1.softbit == 1);
    assert(x1.llr == float.infinity);

    auto x2 = LR01!float(0);
    assert(x2.softbit == -1);
    assert(x2.llr == -float.infinity);

    assert(x1 > x2);
    assert(x2 < LR01!float(1));
    assert(x1 == x1);

    assert((LR01!float(1) & LR01!float(2)) == LR01!float(2));
}


struct LLR01(F, alias mapping = [1, -1])
if(isFloatingPoint!F)
{
    F value;


    static
    LLR01!(F, mapping) fromLR(F lr)
    {
        return LLR01!(F, mapping)(log(lr));
    }


    // compute soft bit
    F softbit() const
    {
        return tanh(value/2);
    }


    typeof(this) toLLR() const
    {
        return this;
    }


    F llr() const
    {
        return value;
    }


    LR01!(F, mapping) toLR() const
    {
        return typeof(return)(exp(value));
    }


    F lr() const
    {
        return this.toLR().value;
    }


    int opCmp(G)(ref const LLR01!(G, mapping) rhs) const
    {
        if(this.value == rhs.value)
            return 0;
        else if(this.value > rhs.value)
            return 1;
        else
            return -1;
    }


    LLR01!(CommonType!(F, G), mapping) opBinary(string op : "+", G)(in LLR01!(G, mapping) rhs) const
    {
        return typeof(return)(this.value + rhs.value);
    }


    LLR01!(CommonType!(F, G), mapping) opBinary(string op : "-", G)(in LLR01!(G, mapping) rhs) const
    {
        return typeof(return)(this.value - rhs.value);
    }


    LLR01!(CommonType!(F, G), mapping) opBinary(string op : "&", G)(in LLR01!(G, mapping) rhs) const
    {
        return this.opBinary!"+"(rhs);
    }


    LLR01!(CommonType!(F, G), mapping) opBinary(string op : "|", G)(in LLR01!(G, mapping) rhs) const
    {
        return this.opBinary!"-"(rhs);
    }
}

unittest
{
    auto x1 = LLR01!float(float.infinity);
    assert(x1.softbit == 1);
    assert(x1.lr == float.infinity);

    auto x2 = LLR01!float(-float.infinity);
    assert(x2.softbit == -1);
    assert(x2.lr == 0);

    assert((LLR01!float(1) & LLR01!float(2)) == LLR01!float(3));
}



struct ProbVector(F, size_t Dim)
if(Dim != 0)
{
    enum F MINIMUM_VALUE = 1E-6;


    F[Dim] values;


    this(in F[] init)
    {
        F sum = 0;
        foreach(i, ref e; values) {
            e = init[i];
        }
        
        this.normalize();
    }


    static
    typeof(this) NO_INFO()
    {
        typeof(return) dst;
        foreach(ref e; dst.values)
            e = 0;

        dst.normalize();
        return dst;
    }


    void normalize()
    {
        F sum = 0;
        foreach(ref e; values) {
            sum += e;
        }

        if(sum < MINIMUM_VALUE) {
            foreach(ref e; values)
                e += MINIMUM_VALUE;
            
            sum += MINIMUM_VALUE * values.length;
        }

        foreach(ref e; values) {
            e /= sum;
        }
    }


    ProbVector!(CommonType!(F, G), Dim) opBinary(string op : "&", G)(in ProbVector!(G, Dim) rhs) const
    {
        typeof(return) dst;
        CommonType!(F, G) sum = 0;
        foreach(i; 0 .. Dim) {
            dst.values[i] = this.values[i] * rhs.values[i];
        }

        dst.normalize();
        return dst;
    }


    ProbVector!(CommonType!(F, G), Dim) opBinary(string op : "|", G)(in ProbVector!(G, Dim) rhs) const
    {
        typeof(return) dst;
        CommonType!(F, G) sum = 0;
        foreach(i; 0 .. Dim) {
            dst.values[i] = this.values[i] / rhs.values[i];
        }


        dst.normalize();
        return dst;
    }


    void opOpAssign(string op, G)(in ProbVector!(G, Dim) rhs)
    {
        auto vs = this.opBinary!op(rhs);
        this.values[] = vs.values[];
    }


    G expect(G)(in G[] arr) const
    {
        G sum = 0;
        foreach(i; 0 .. Dim)
            sum += arr[i] * values[i];

        return sum;
    }
}


unittest
{
    import std.stdio;
    auto x1 = ProbVector!(float, 4).NO_INFO;
    assert(x1.expect([1.0, 2, 3, 4]).isClose((1 + 2 + 3 + 4) / 4.0));

    x1 &= ProbVector!(float, 4)([0, 0, 1, 0]);
    // writeln(x1.expect([1.0, 2, 3, 4]));
    assert(x1.expect([1.0, 2, 3, 4]).isClose(3));

    x1 &= ProbVector!(float, 4)([1, 0, 0, 0]);
    assert(x1 == typeof(x1).NO_INFO);
}


// unittest
// {
//     auto x1 = LLR01!float(float.infinity);
//     assert(x1.softbit == 1);
//     assert(x1.lr == float.infinity);

//     auto x2 = LLR01!float(-float.infinity);
//     assert(x2.softbit == -1);
//     assert(x2.lr == 0);

//     assert((LLR01!float(1) & LLR01!float(2)) == LLR01!float(3));
// }


/**
 事後平均とその微分値と事後分散を計算します
*/
Tuple!(F, "value", F, "drv", F, "var") softDecision(alias arrP_, alias arrR_, F)(F x, F sigma2)
if(isArray!(typeof(arrP_)) && isArray!(typeof(arrR_)) && isFloatingPoint!F && arrP_.length >= 2 && arrP_.length == arrR_.length)
{
    import std.math : isNaN, abs;
    import std.algorithm : min;
    import dffdd.math.math : fast_exp;

    enum F[] arrP = arrP_;
    enum F[] arrR = arrR_;

    F xr2_min = F.infinity;
    static foreach(i; 0 .. arrP.length) {
        if(arrP[i] != 0)
            xr2_min = min((arrR[i] - x)^^2, xr2_min);
    }

    immutable expCoef = -0.5 / sigma2;

    F p_e = F(0),
    p_r_e = F(0),
    p_r2_e = F(0);

    F[arrP.length] _pes;
    static foreach(i; 0 .. arrP.length) {{
        if(arrP[i] != 0) {
            immutable xr2 = (arrR[i] - x)^^2;

            F expvalue = fast_exp!F(expCoef * (xr2 - xr2_min));
            if(expvalue.isNaN) expvalue = 1;

            _pes[i] = arrP[i] * expvalue;
            p_e += _pes[i];
            p_r_e += _pes[i] * arrR[i];
            p_r2_e += _pes[i] * arrR[i]^^2;
        }
    }}

    immutable drv = -2 * expCoef * (p_r2_e * p_e - p_r_e^^2) / p_e^^2;
    immutable mean = p_r_e / p_e;

    return typeof(return)(
        mean,
        (abs(sigma2) > 1e-30 && !drv.isNaN) ? drv : 0,
        (p_r2_e / p_e) - mean^^2
    );
}

unittest
{
    assert(softDecision!([0.5, 0, 0.5], [-1, 0, 1])(0.001, (0.1/sqrt(0.8))^^2 ).value.isClose(0.0798298, 1e-4));
    assert(softDecision!([0.5, 0, 0.5], [-1, 0, 1])(0.001, (0.1/sqrt(0.8))^^2 ).drv.isClose(79.4902, 1e-4));
    assert(softDecision!([0.5, 0, 0.5], [-1, 0, 1])(0.001, (0.1/sqrt(0.8))^^2 ).var.isClose(1 - 0.0798298^^2, 1e-4));   
}


alias softDecision2PAM_BPSK(F) = softDecision!([0.5, 0.5], [-1, 1], F);
alias softDecision2PAM_QPSK(F) = softDecision!([0.5, 0.5], [-SQRT1_2, SQRT1_2], F);
alias softDecision4PAM_16QAM(F) = softDecision!([0.25, 0.25, 0.25, 0.25], [-3/sqrt(10.0L), -1/sqrt(10.0L), 1/sqrt(10.0L), 3/sqrt(10.0L)], F);
alias softDecision8PAM_64QAM(F) = softDecision!([0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125], [-7/sqrt(42.0L), -5/sqrt(42.0L), -3/sqrt(42.0L), -1/sqrt(42.0L), 1/sqrt(42.0L), 3/sqrt(42.0L), 5/sqrt(42.0L), 7/sqrt(42.0L)], F);
