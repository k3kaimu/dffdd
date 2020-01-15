module dffdd.math.math;

import std.complex : Complex;
import std.math;
import std.traits;


/**
compute x^^y
*/
F libm_pow(F)(F x, F y) nothrow @nogc
if(isFloatingPoint!F)
{
    import core.stdc.math : powf, pow, powl;

    static if(is(F == float))
    {
        return powf(x, y);
    }
    else static if(is(F == double))
    {
        return pow(x, y);
    }
    else
    {
        return powl(x, y);
    }
}


/**
compute log(x)
*/
F libm_log(F)(F x, F y) nothrow @nogc
if(isFloatingPoint!F)
{
    import core.stdc.math : logf, log, logl;

    static if(is(F == float))
    {
        return logf(x);
    }
    else static if(is(F == double))
    {
        return log(x);
    }
    else
    {
        return logl(x);
    }
}


/**
compute sin(x)
*/
F libm_sin(F)(F x) nothrow @nogc
if(isFloatingPoint!F)
{
    import core.stdc.math : sinf, sin, sinl;

    static if(is(F == float))
    {
        return sinf(x);
    }
    else static if(is(F == double))
    {
        return sin(x);
    }
    else
    {
        return sinl(x);
    }
}


/**
compute cos(x)
*/
F libm_cos(F)(F x) nothrow @nogc
if(isFloatingPoint!F)
{
    import core.stdc.math : cosf, cos, cosl;

    static if(is(F == float))
    {
        return cosf(x);
    }
    else static if(is(F == double))
    {
        return cos(x);
    }
    else
    {
        return cosl(x);
    }
}



version(CRuntime_Glibc)
extern(C)
{
    void sincosf(float x, float *sin, float *cos) nothrow @nogc;
    void sincos(double x, double *sin, double *cos) nothrow @nogc;
    void sincosl(real x, real *sin, real *cos) nothrow @nogc;
}


/**
compute exp(ix)
*/
Complex!F libm_expi(F)(F x) nothrow @nogc
if(isFloatingPoint!F)
{
    F re, im;

    static if(is(F == float) && is(typeof(&sincosf)))
    {
        sincosf(x, &im, &re);
    }
    else static if(is(F == double) && is(typeof(&sincos)))
    {
        sincos(x, &im, &re);
    }
    else static if(is(F == real) && is(typeof(&sincosl)))
    {
        sincosl(x, &im, &re);
    }
    else
    {
        re = libm_cos(x);
        im = libm_sin(x);
    }

    return Complex!F(re, im);
}


/**
compute sqrt(x)
*/
F libm_sqrt(F)(F x) nothrow @nogc
if(isFloatingPoint!F)
{
    import core.stdc.math : sqrtf, sqrt, sqrtl;

    static if(is(F == float))
    {
        return sqrtf(x);
    }
    else static if(is(F == double))
    {
        return sqrt(x);
    }
    else
    {
        return sqrtl(x);
    }
}


/**
compute cbrt(x)
*/
F libm_cbrt(F)(F x) nothrow @nogc
if(isFloatingPoint!F)
{
    import core.stdc.math : cbrtf, cbrt, cbrtl;

    static if(is(F == float))
    {
        return cbrtf(x);
    }
    else static if(is(F == double))
    {
        return cbrt(x);
    }
    else
    {
        return cbrtl(x);
    }
}


/**
compute log(x)
*/
F libm_log(F)(F x) nothrow @nogc
if(isFloatingPoint!F)
{
    import core.stdc.math : logf, log, logl;

    static if(is(F == float))
    {
        return logf(x);
    }
    else static if(is(F == double))
    {
        return log(x);
    }
    else
    {
        return logl(x);
    }
}


version(LDC)
{
    pragma(LDC_intrinsic, "llvm.pow.f#")
    F llvm_pow(F)(F x, F y) nothrow @nogc;
    pragma(LDC_intrinsic, "llvm.powi.f#")
    F llvm_powi(F)(F x, int y) nothrow @nogc;
    pragma(LDC_intrinsic, "llvm.sqrt.f#")
    F llvm_sqrt(F)(F x) nothrow @nogc;
    pragma(LDC_intrinsic, "llvm.sin.f#")
    F llvm_sin(F)(F x) nothrow @nogc;
    pragma(LDC_intrinsic, "llvm.cos.f#")
    F llvm_cos(F)(F x) nothrow @nogc;
    pragma(LDC_intrinsic, "llvm.log.f#")
    F llvm_log(F)(F x) nothrow @nogc;

    
    Complex!F llvm_expi(F)(F x) nothrow @nogc
    if(isFloatingPoint!F)
    {
        F re = llvm_cos(x);
        F im = llvm_sin(x);
        return Complex!F(re, im);
    }


    alias fast_pow = llvm_pow;
    alias fast_powi = llvm_powi;
    alias fast_sqrt = llvm_sqrt;
    alias fast_cbrt = libm_cbrt;
    alias fast_sin = llvm_sin;
    alias fast_cos = llvm_cos;
    alias fast_expi = llvm_expi;
    alias fast_log = llvm_log;
}
else
{
    alias fast_pow = libm_pow;
    alias fast_sqrt = libm_sqrt;
    alias fast_cbrt = libm_cbrt;
    alias fast_sin = libm_sin;
    alias fast_cos = libm_cos;
    alias fast_expi = libm_expi;
    alias fast_log = libm_log;

    F fast_powi(F)(F x, int y) nothrow @nogc
    {
        return x^^y;
    }
}


/**
compute |x|
*/
auto fast_abs(C)(C x)
{
    static if(is(typeof(C.init.re)))
    {
        // for complex number
        alias F = typeof(C.init.re);
        F p = x.re*x.re + x.im*x.im;
        return fast_sqrt!F(p);
    }
    else
    {
        // for real number
        return cast(C)std.math.abs(x);
    }
}


/**
0次の変形ベッセル関数I_0(x)を計算します．
*/
F besselI0(F)(F x) pure nothrow @safe @nogc
if(isFloatingPoint!F)
{
    static if(is(F == const) || is(F == immutable))
        return besselI0!(Unqual!F)(x);
    else
    {
        if(x.isNaN) return F.nan;
        if(x.isInfinity) return x > 0 ? F.infinity : F.nan;

        x *= 0.5;

        F dst = 0;
        F r = 1;
        size_t cnt = 0;
        while(abs(r) > abs(F.epsilon * dst)){
            dst += r;
            ++cnt;

            immutable s = x / cnt;
            r *= s * s;
        }

        return dst;
    }
}

unittest
{
    // compare with numpy's result
    assert(approxEqual(besselI0(0.0L),      1));
    assert(approxEqual(besselI0(1.0L),      1.2660658777520082));
    assert(approxEqual(besselI0(100.0L),    1.0737517071310738e+42));

    assert(approxEqual(besselI0(0.1)-1,     0.002501562934095647));
    assert(approxEqual(besselI0(0.01)-1,    2.5000156250509775e-05));
}


/**
sin(pi * x) / (pi * x)を計算します
*/
F sinc(F)(F x) pure nothrow @safe @nogc
if(isFloatingPoint!F)
{
    if(abs(x) < 1E-3){
        immutable px2 = (PI * x)^^2;
        F[5] coefs = [
            1,
            -1/(2.0L*3),
            1/(2.0L*3*4*5),
            -1/(2.0L*3*4*5*6*7),
            1/(2.0L*3*4*5*6*7*8*9),
        ];

        return poly(px2, coefs);
    }else{
        immutable px = PI * x;
        return sin(px) / px;
    }
}

unittest
{
    // compare with numpy's result
    assert(approxEqual(1-sinc(0.001),   1.6449332550516615e-06));
    assert(approxEqual(1-sinc(0.0001),  1.644934066735715e-08));
    assert(approxEqual(1-sinc(0.0002),  6.579736144818327e-08));
    assert(approxEqual(sinc(0.1),       0.983631643083466));
}


/**
超幾何級数 2F1(a, b; c; z) を計算します
See also: https://ja.wikipedia.org/wiki/%E8%B6%85%E5%B9%BE%E4%BD%95%E7%B4%9A%E6%95%B0
*/
F hypgeom21(F = real)(F a, F b, F c, F z)
{
    if(a.isNaN || b.isNaN || c.isNaN || z.isNaN)
        return F.nan;
    if((a > 0 && a.isInfinity) || (b > 0 && b.isInfinity) || (z > 0 && z.isInfinity))
        return F.infinity;
    if((a < 0 && a.isInfinity) || (b < 0 && b.isInfinity) || (z < 0 && z.isInfinity))
        return F.nan;

    F sum = 0;
    F next = 1;
    size_t cnt = 0;
    while(abs(next) > abs(F.epsilon * sum) || sum.isInfinity || sum.isNaN){
        sum += next;

        ++cnt;
        next = next * (a + (cnt-1)) * (b + (cnt-1)) / cnt / (c + (cnt-1)) * z;
    }

    return sum;
}

unittest 
{
    real mylog(real x){
        x -= 1;
        return x * hypgeom21(1, 1, 2, -x);
    }

    assert(approxEqual(mylog(1.01), log(1.01)));
    assert(approxEqual(mylog(1.1), log(1.1)));
}


/**
\Lambda_0^{*}(x) = int [0 -> x] I0(t) dt を計算します
See also: http://web.eah-jena.de/~rsh/Forschung/Stoer/besint.pdf page: 125
*/
F intBesselI0(F = real)(F x)
{
    if(x.isNaN) return F.nan;
    if(x.isInfinity) return x > 0 ? F.infinity : F.nan;

    immutable ox = x;
    x *= 0.5;

    F dst = 0;
    F r = 1;
    size_t cnt = 0;
    while(abs(r) > abs(F.epsilon * dst)){
        dst += r / (2*cnt + 1);
        ++cnt;

        immutable s = x / cnt;
        r *= s * s;
    }

    return ox * dst;
}

unittest 
{
    // see also: http://web.eah-jena.de/~rsh/Forschung/Stoer/besint.pdf p. 123
    real[] alphas = [
        1.0000000000000000E+00,
        -8.3333333333333329E-02,
        3.1250000000000002E-03,
        -6.2003968253968251E-05,
        7.5352044753086416E-07,
        -6.1651672979797980E-09,
        3.6226944592830012E-11,
        -1.6018716996829596E-13,
        5.5211570531352012E-16,
        -1.5246859958300586E-18,
        3.4486945143775135E-21,
        -6.5058017249306315E-24,
        1.0391211088430869E-26,
        -1.4232975959389203E-29,
        1.6902284962328838E-32,
        -1.7568683294176930E-35,
        1.6117104111016952E-38,
        -1.3145438350557573E-41,
        9.5948102742224525E-45,
        -6.3038564554696844E-48,
        3.7477195390749646E-51,
    ];

    foreach(ref e; alphas) e = abs(e);

    real approxLambda0ast(real x)
    {
        return x * poly(x^^2, alphas);
    }

    assert(approxEqual(approxLambda0ast(0.1), intBesselI0(0.1)));
}