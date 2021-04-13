module dffdd.blockdiagram.filter;


import std.range;
import std.traits;
import std.json;

import carbon.math : nextPowOf2;

import dffdd.utils.json;
import dffdd.utils.unit;


struct FIRFilterConverter(C)
{
    import dffdd.utils.fft;
    import std.complex;
    import std.algorithm : min;

    alias F = typeof(C.init.re);

    alias InputElementType = C;
    alias OutputElementType = C;


    this(in C[] coefs)
    {
        _coefs = coefs.dup;
        _inputs.length = coefs.length;

        foreach(ref e; _inputs)
            e = C(0);

        immutable size = nextPowOf2(_coefs.length) * 4;

        _fftw = makeFFTWObject!Complex(size);

        // 係数をFFT
        _fftw.inputs!F[0 .. _coefs.length] = _coefs[];
        _fftw.inputs!F[_coefs.length .. size] = C(0);
        _fftw.fft!F();
        _freqCoefs = _fftw.outputs!F.dup;
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        foreach(i; 0 .. _inputs.length-1)
            _inputs[i] = _inputs[i+1];

        _inputs[$-1] = input;

        output = C(0);
        foreach(i; 0 .. _coefs.length)
            output += _coefs[$-1-i] * _inputs[i];
    }


    void opCall(const(InputElementType)[] input, ref OutputElementType[] output)
    {
        if(output.length != input.length)
            output.length = input.length;

        immutable clen = _coefs.length;
        immutable size = _freqCoefs.length;

        auto dst = output;

        while(input.length != 0){
            immutable ss = min(size - (clen-1), input.length);

            _fftw.inputs!F[0 .. clen-1] = _inputs[1 .. $];
            _fftw.inputs!F[clen-1 .. clen-1 + ss] = input[0 .. ss];
            _fftw.inputs!F[clen-1 + ss .. $] = C(0);
            _fftw.fft!F();

            // 係数と入力信号の積
            _fftw.inputs!F[] = _freqCoefs[] * _fftw.outputs!F[];
            _fftw.ifft!F();

            dst[0 .. ss] = _fftw.outputs!F[clen-1 .. clen-1 + ss];
            
            // 今回使用した入力信号の量がタップ数より多い場合
            if(ss >= clen){
                // 後ろからclen分だけコピー
                _inputs[] = input[ss - clen .. ss];
            }else{
                // ss分だけ_inputを進める
                foreach(i; 0 .. clen-ss)
                    _inputs[i] = _inputs[i+ss];

                _inputs[$-ss .. $] = input[0 .. ss];
            }

            input = input[ss .. $];
            dst = dst[ss .. $];
        }
    }


    FIRFilterConverter dup() const @property
    {
        typeof(return) dst;
        dst._coefs = this._coefs.dup;
        dst._inputs = this._inputs.dup;
        dst._fftw = makeFFTWObject!Complex(this._freqCoefs.length);
        dst._freqCoefs = this._freqCoefs.dup;

        return dst;
    }


    const(C)[] coefficients() const @property
    {
        return _coefs;
    }


    JSONValue dumpInfoToJSON() const
    {
        return JSONValue([
            "coefs": DefaultJSONEnv.toJSONValue(_coefs)
        ]);
    }


  private:
    immutable(C)[] _coefs;
    C[] _inputs;
    FFTWObject!Complex _fftw;
    immutable(C)[] _freqCoefs;
}

unittest
{
    import dffdd.blockdiagram.utils;
    import std.algorithm;
    import std.complex;
    import std.math;
    import std.range;
    import std.stdio;

    alias C = Complex!float;

    auto sig = [C(0), C(1), C(2), C(3)].connectTo!(FIRFilterConverter!C)([C(0), C(1)]);
    assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(sig, [C(0), C(0), C(1), C(2)]));
}



struct FIRFilter
{
    static
    auto makeBlock(R, C)(R r, C[] coefs)
    if(isInputRange!R)
    {
        import dffdd.blockdiagram.utils;
        alias E = Unqual!(ElementType!R);
        return r.connectTo!(FIRFilterConverter!E)(coefs);
    }
}

unittest
{
    import std.algorithm;
    import std.range;
    import std.stdio;
    import std.complex;

    alias C = Complex!float;

    auto sig = FIRFilter.makeBlock([C(0), C(1), C(2), C(3)], [C(0), C(1)]);
    assert(equal(sig, [C(0), C(0), C(1), C(2)]));
}


struct IIRFilter
{
    static
    auto makeBlock(R, C)(R r, C coef1)
    if(isInputRange!R)
    {
        return Simple1TapIIRFilter!(R, C)(r, coef1);
    }


    static struct Simple1TapIIRFilter(R, C)
    {
        private alias E = Unqual!(ElementType!R);
        private alias F = typeof(C.init * E.init);

        this(R)(R r, C coef1)
        {
            _r = r;
            _coef1 = coef1;
            _front = cast(F)0;
            popFront();
        }


        F front() const @property { return _front; }
        bool empty() const @property { return _empty; }


        void popFront()
        {
            if(_r.empty){
                _empty = true;
                return;
            }

            _front = _r.front + _coef1 * _front;
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
        C _coef1;
        F _front;
        bool _empty;
    }
}
