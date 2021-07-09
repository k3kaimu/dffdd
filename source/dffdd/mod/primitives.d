module dffdd.mod.primitives;

import std.numeric : gcd;
import std.typecons;



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






struct P0P1Calculator(C, F)
{
    this(Mod)(Mod mod)
    in(mod.symOutputLength == 1)
    {
        _nbits = mod.symInputLength;
        _p0 = new F[_nbits];
        _p1 = new F[_nbits];

        C[] sym = new C[1];
        Bit[] bits = new Bit[_nbits];
        foreach(n; 0 .. 1 << _nbits) {
            foreach(k; 0 .. _nbits)
                bits[k] = (n & (1 << k)) ? 1 : 0;

            mod.modulate(bits, sym);
            _points ~= sym[0];
        }
    }


    ref F[] computeP0P1(in C[] inputs, return ref F[] p0p1, F N0)
    {
        p0p1.length = inputs.length * _nbits;

        foreach(i, x; inputs) {
            _p0[] = 0;
            _p1[] = 0;

            foreach(j, y; _points) {
                auto p = exp(-(x - y).sqAbs/ N0);

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
        
        auto bits = new Bit[2^^_nbits * _nbits];
        foreach(n; 0 .. 2^^_nbits) {
            foreach(j; 0 .. _nbits) {
                bits ~= Bit(n & (1 << j));
            }
        }

        _symbols = mod.modulate(bits, _symbols);
    }


    ref C[] opCall(in R[] p0p1s, return ref C[] outputs)
    in(p0p1s.length % _nbits == 0)
    do {
        outputs.length = p0p1s.length / _nbits;

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

