module dffdd.blockdiagram.fft;

import std.typecons;

import std.complex;
import carbon.complex;

import dffdd.blockdiagram.utils;
import rx;

import dffdd.utils.fft;


struct FastFourierTransformer(C, Flag!"isForward" isForward)
if(isComplex!C)
{
    alias InputElementType = C;
    alias OutputElementType = C;

    alias F = typeof(C.init.re);


    this(size_t n)
    {
        _fftw = makeFFTWObject!C(n);
        _size = n;
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    in{
        assert(input.length == _size);
    }
    body{
        _fftw.inputs!F[] = input[];
        
        static if(isForward)
            _fftw.fft!F();
        else
            _fftw.ifft!F();
        
        output.length = _size;
        output[] = _fftw.outputs!F[];
    }


    FastFourierTransformer dup() const
    {
        return FastFourierTransformer(_size);
    }


  private:
    FFTWObject!C _fftw;
    size_t _size;
}


auto makeFastFourierTransformer(C)(size_t n)
{
    return FastFourierTransformer!(C, Yes.isForward)(n);
}


auto makeInverseFastFourierTransformer(C)(size_t n)
{
    return FastFourierTransformer!(C, No.isForward)(n);
}


unittest 
{
    import std.random;
    import std.stdio : writeln;
    import std.range;
    import std.algorithm;
    import std.math;

    alias C = Complex!float;

    C[] rnds;
    foreach(i; 0 .. 1024)
        rnds ~= C(uniform01, uniform01);

    auto vals = rnds.chunks(32)
        .connectTo(makeFastFourierTransformer!C(32))
        .connectTo(makeInverseFastFourierTransformer!C(32)).joiner();

    assert(equal!((a, b) => approxEqual(a.re, b.re) && approxEqual(a.im, b.im))(rnds, vals));
}
