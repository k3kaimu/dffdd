module dffdd.mod.primitives;

import std.numeric : gcd;



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
    body{
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
    body{
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
    import std.complex;
    import dffdd.mod.ofdm;
    import dffdd.mod.qpsk;

    // 1次変調: QPSK
    // FFTサイズ: 8
    // サイクリックプレフィックス: 3
    // 使用サブキャリア数: 4
    // アップサンプリング率: 2
    auto mod = chainedMod(new QPSK(), new OFDM!(Complex!float)(8, 3, 4, 2));

    assert(mod.symInputLength == 2 * 4);
    assert(mod.symOutputLength == (8 + 3) * 2);

    ubyte[] inps = [1, 1, 0, 0, 1, 0, 0, 1];
    Complex!float[] res;

    // 変調
    mod.modulate(inps, res);

    // CPのチェック
    foreach(i; 0 .. 3 * 2)
        assert(res[i] == res[$ - 3*2 + i]);

    ubyte[] decoded;

    // 復調
    mod.demodulate(res, decoded);

    // 復調結果のチェック
    assert(decoded == inps);
}
