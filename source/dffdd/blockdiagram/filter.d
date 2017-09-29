module dffdd.blockdiagram.filter;


import std.range;
import std.traits;

import carbon.math : nextPowOf2;


struct FIRFilter
{
    static
    auto makeBlock(R, C)(R r, C[] coefs)
    if(isInputRange!R)
    {
        return SimpleFIRFilter!(R, C)(r, coefs);
    }


    static
    auto makeBlock(R, C)(R r, C coef1)
    if(isInputRange!R && is(typeof((R r, C c){ auto e = r.front * c; })))
    {
        return Simple1TapFIRFilter!(R, C)(r, coef1);
    }



    static struct SimpleFIRFilter(R, C)
    {
        private alias E = Unqual!(ElementType!R);
        private alias F = typeof(C.init * E.init);


        this(R r, C[] coefs)
        {
            _r = r;
            _coefs = coefs.dup;
            _regs = (new E[coefs.length]).cycle();

            foreach(i; 0 .. coefs.length)
                _regs[i] = cast(F)0;

            this.popFront();
        }


        F front() const @property { return _front; }
        bool empty() const @property { return _empty; }

        void popFront()
        {
            _regs.popFront();

            if(_r.empty){
                _empty = true;
                return;
            }
            _regs[_coefs.length - 1] = _r.front;
            _r.popFront();
            _front = cast(F)0;
            foreach(i; 0 .. _coefs.length)
                _front += _coefs[$-1-i] * _regs[i];
        }


      static if(isForwardRange!R)
      {
        typeof(this) save() @property
        {
            typeof(return) dst = this;

            dst._r = this._r.save;
            dst._regs = (new E[this._coefs.length]).cycle;
            foreach(i; 0 .. this._coefs.length)
                dst._regs[i] = this._regs[i];

            return dst;
        }
      }


      private:
        R _r;
        immutable(C)[] _coefs;
        typeof(cycle(new E[1])) _regs;
        F _front;
        bool _empty;
    }


    static struct Simple1TapFIRFilter(R, C)
    {
        private alias E = Unqual!(ElementType!R);
        private alias F = typeof(C.init * E.init);


        this(R)(R r, C coef1)
        {
            _r = r;
            _coef1 = coef1;
            _x0 = cast(E)0;
            _x1 = cast(E)0;
            popFront();
        }


        F front() const @property { return _front; }
        bool empty() const @property { return _empty; }


        void popFront()
        {
            _x1 = _x0;

            if(_r.empty){
                _empty = true;
                return;
            }

            _x0 = _r.front;
            _r.popFront();

            _front = _x0 + _x1*_coef1;
        }


      static if(isForwardRange!R)
      {
        typeof(this) save() @property
        {
            typeof(return) dst = this;

            dst._r = this._r.save;
        }
      }


      private:
        R _r;
        C _coef1;
        E _x0, _x1;
        F _front;
        bool _empty;
    }
}

unittest
{
    import std.algorithm;
    import std.range;
    import std.stdio;

    auto sig = FIRFilter.makeBlock([0, 1, 2, 3], [0, 1]);
    assert(equal(sig, [0, 0, 1, 2]));
}


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
        _coefs = coefs;
        _inputs.length = coefs.length;

        foreach(ref e; _inputs)
            e = C(0);

        immutable size = nextPowOf2(_coefs.length) * 4;

        _fftw = makeFFTWObject!Complex(size);
        _tempbuffer = new C[size];
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
        import std.stdio : writeln;

        immutable clen = _coefs.length;
        immutable size = _tempbuffer.length;

        output.length = input.length;

        // 係数をFFT
        _fftw.inputs!F[0 .. clen] = _coefs[];
        _fftw.inputs!F[clen .. size] = C(0);
        _fftw.fft!F();
        _tempbuffer[] = _fftw.outputs!F[];

        auto dst = output;

        while(input.length != 0){
            immutable ss = min(size - (clen-1), input.length);

            // writeln(__LINE__);
            _fftw.inputs!F[0 .. clen-1] = _inputs[1 .. $];
            // writeln(__LINE__);
            _fftw.inputs!F[clen-1 .. clen-1 + ss] = input[0 .. ss];
            // writeln(__LINE__);
            _fftw.inputs!F[clen-1 + ss .. $] = C(0);
            _fftw.fft!F();

            // 係数と入力信号の積
            _fftw.inputs!F[] = _tempbuffer[] * _fftw.outputs!F[];
            _fftw.ifft!F();

            // writeln(__LINE__);
            dst[0 .. ss] = _fftw.outputs!F[clen-1 .. clen-1 + ss];
            
            // 今回使用した入力信号の量がタップ数より多い場合
            if(ss >= clen){
                // 後ろからclen分だけコピー
                // writeln(__LINE__);
                _inputs[] = input[ss - clen .. ss];
            }else{
                // ss分だけ_inputを進める
                foreach(i; 0 .. clen-ss)
                    _inputs[i] = _inputs[i+ss];
                // writeln(__LINE__);
                _inputs[$-ss .. $] = input[0 .. ss];
            }

            // writeln(__LINE__);
            input = input[ss .. $];
            dst = dst[ss .. $];
        }
    }


    FIRFilterConverter dup() const @property
    {
        typeof(return) dst;
        dst._coefs = this._coefs;
        dst._inputs = this._inputs.dup;
        dst._fftw = makeFFTWObject!Complex(this._tempbuffer.length);
        dst._tempbuffer = this._tempbuffer.dup;

        return dst;
    }


  private:
    const(C)[] _coefs;
    C[] _inputs;
    FFTWObject!Complex _fftw;
    C[] _tempbuffer;
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

unittest 
{
    import dffdd.blockdiagram.utils;
    import std.algorithm;
    import std.complex;
    import std.random;
    import std.range;

    alias C = Complex!float;

    // static assert(0, "You should write test code!!!!");
    static
    void testImpl(size_t ntap, size_t n)
    {
        import std.stdio : writefln, writeln;
        import std.math;

        C[] taps = new C[ntap];
        foreach(i; 0 .. ntap)
            taps[i] = C(uniform01(), uniform01());

        auto conv1 = FIRFilterConverter!C(taps.dup);
        auto conv2 = FIRFilterConverter!C(taps.dup);

        C[] buffers = new C[n * 4];
        foreach(i; 0 .. n * 4)
            buffers[i] = C(uniform01(), uniform01());

        auto r0 = FIRFilter.makeBlock(buffers, taps.dup).array();
        auto r1 = buffers.chunks(n).connectTo!(FIRFilterConverter!C)(taps.dup).joiner;
        auto r2 = buffers.connectTo!(FIRFilterConverter!C)(taps.dup);

        assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r1));
        assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(r0, r2));
    }

    testImpl(16, 64);
    testImpl(15, 64);
    testImpl(12, 64);
    testImpl(16, 16);
    testImpl(16, 17);
    testImpl(16, 1);
    testImpl(16, 2);
    testImpl(16, 3);
    testImpl(5, 3);
    testImpl(5, 9);
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
