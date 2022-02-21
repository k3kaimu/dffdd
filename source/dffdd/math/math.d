module dffdd.math.math;

// import std.complex;
import std.math;
import std.traits;

import dffdd.math.complex;

/**
compute x^^y
*/
F libm_pow(F)(F x, F y) /*pure*/ nothrow @nogc
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
compute sqrt(x)
*/
F libm_sqrt(F)(F x) /*pure*/ nothrow @nogc
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
F libm_cbrt(F)(F x) /*pure*/ nothrow @nogc
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
compute exp(x)
*/
F libm_exp(F)(F x) /*pure*/ nothrow @nogc
if(isFloatingPoint!F)
{
    import core.stdc.math : expf, exp, expl;

    static if(is(F == float))
    {
        return expf(x);
    }
    else static if(is(F == double))
    {
        return exp(x);
    }
    else
    {
        return expl(x);
    }
}


version(LDC)
{
    pragma(LDC_intrinsic, "llvm.pow.f#")
    F llvm_pow(F)(F x, F y) pure nothrow @nogc;
    pragma(LDC_intrinsic, "llvm.powi.f#")
    F llvm_powi(F)(F x, int y) pure nothrow @nogc;
    pragma(LDC_intrinsic, "llvm.sqrt.f#")
    F llvm_sqrt(F)(F x) pure nothrow @nogc;
    pragma(LDC_intrinsic, "llvm.exp.f#")
    F llvm_exp(F)(F x) pure nothrow @nogc;


    alias fast_pow = llvm_pow;
    alias fast_powi = llvm_powi;
    alias fast_sqrt = llvm_sqrt;
    alias fast_cbrt = libm_cbrt;
    alias fast_exp = llvm_exp;
}
else
{
    alias fast_pow = libm_pow;
    alias fast_sqrt = libm_sqrt;
    alias fast_cbrt = libm_cbrt;
    alias fast_exp = libm_exp;

    F fast_powi(F)(F x, int y) pure nothrow @nogc
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
    assert(isClose(besselI0(0.0L),      1));
    assert(isClose(besselI0(1.0L),      1.2660658777520082));
    assert(isClose(besselI0(100.0L),    1.0737517071310738e+42));

    assert(isClose(besselI0(0.1)-1,     0.002501562934095647));
    assert(isClose(besselI0(0.01)-1,    2.5000156250509775e-05));
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
    assert(isClose(1-sinc(0.001),   1.6449332550516615e-06, 1e-4));
    assert(isClose(1-sinc(0.0001),  1.644934066735715e-08, 1e-4));
    assert(isClose(1-sinc(0.0002),  6.579736144818327e-08, 1e-4));
    assert(isClose(sinc(0.1),       0.983631643083466, 1e-4));
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

    assert(isClose(mylog(1.01), log(1.01)));
    assert(isClose(mylog(1.1), log(1.1)));
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

    assert(isClose(approxLambda0ast(0.1), intBesselI0(0.1)));
}




C fact(C = ulong)(uint n)
{
    C dst = 1;
    foreach(i; 2 .. n + 1)
        dst *= i;

    return dst;
}

unittest
{
    assert(fact(0) == 1);
    assert(fact(1) == 1);
    assert(fact(2) == 2);
    assert(fact(3) == 6);
    assert(fact(4) == 24);
    assert(fact(5) == 120);
}

C binom(C)(C x, size_t k) if(isComplex!C || isFloatingPoint!C || isIntegral!C)
{
    C num = 1;
    C den = 1;
    foreach(i; 0 .. k)
    {
        num *= x - i;
        den *= k - i;
    }

    return num / den;
}

unittest
{
    assert(binom(5, 0) == 1);
    assert(binom(5, 1) == 5);
    assert(binom(5, 2) == 10);
    assert(binom(5, 3) == 10);
    assert(binom(5, 4) == 5);
    assert(binom(5, 5) == 1);
}

/**
一般ラゲール多項式（Generalized Laguerre Polynomial） L_n^a(x)を評価します．
*/
C laguerre(C)(C x, int n, C a = 0) if(isComplex!C || isFloatingPoint!C)
{
    C sum = 0;
    foreach(i; 0 .. n + 1)
        sum += (-1) ^^ i * binom(n + a, n - i) * x ^^ i / fact!C(i);

    return sum;
}

unittest
{
    import std.math : isClose;

    // L_0^a(x) = 1
    auto l0a = delegate(real a, real x) { return 1.0L; };

    // L_1^a(x) = -x + a + 1
    auto l1a = delegate(real a, real x) { return -x + a + 1; };

    // L_2^a(x) = x^2/2 - (a+2)x + (a+2)(a+1)/2
    auto l2a = delegate(real a, real x) {
        return x ^^ 2 / 2 - (a + 2) * x + (a + 2) * (a + 1) / 2;
    };

    // L_3^a(x) = -x^3/6 + (a+3)x^2/2 - (a+2)(a+3)x/2 + (a+1)(a+2)(a+3)/6
    auto l3a = delegate(real a, real x) {
        return -x ^^ 3 / 6 + (a + 3) * x ^^ 2 / 2 - (a + 2) * (a + 3) * x / 2 + (
                a + 1) * (a + 2) * (a + 3) / 6;
    };

    auto testXs = [0, 0.1, 0.5, 1, 2, 5, 10];
    auto testAs = [0, 1, 2, 3, 4, 5];

    foreach(x; testXs)
        foreach(a; testAs)
        {
            foreach(int i, f; [l0a, l1a, l2a, l3a])
                assert(laguerre(x, i, a).isClose(f(a, x)));
        }
}

/**
一般ラゲール多項式に基づく正規直交基底を計算します
*/
C laguerreOBF(C)(C x, int m, int n) if(isComplex!C || isFloatingPoint!C)
{
    alias F = typeof(abs(x));

    if(m >= n)
    {
        immutable F normCoef = sqrt(fact!F(n) / fact!F(m));
        return (-1) ^^ n * normCoef * x ^^ (m - n) * laguerre(x.sqAbs, n, m - n);
    }
    else
    {
        immutable F normCoef = sqrt(fact!F(m) / fact!F(n));
        return (-1) ^^ m * normCoef * x.conj ^^ (n - m) * laguerre(x.sqAbs, m, n - m);
    }
}

unittest
{
    //
    auto ps = [0, 1, 2, 3, 4, 5];
    auto xs = [stdComplex(0, 0), stdComplex(0, 1), stdComplex(1, 0), stdComplex(1, 1),
        stdComplex(0.5, 0.5), stdComplex(-0.5, -0.5), stdComplex(0.5, -0.5), stdComplex(-0.5, 0.5)];

    foreach(p; ps)
    {
        foreach(x; xs)
        {
            assert(laguerreOBF(x, p, 0).approxEqualC(1 / sqrt(fact!real(p)) * x ^^ p));
            assert(laguerreOBF(x, 0, p).approxEqualC(laguerreOBF(x, p, 0).conj));
            assert(laguerreOBF(x, p, p).approxEqualC((-1) ^^ p * laguerre(x.sqAbs, p)));
            assert(laguerreOBF(x, p + 1, p).approxEqualC((-1) ^^ p / sqrt(p + 1.0L) * x * laguerre(x.sqAbs,
                    p, 1)));
            assert(laguerreOBF(x, p, p + 1).approxEqualC(laguerreOBF(x, p + 1, p).conj));
        }
    }
}


bool approxEqualC(C1, C2)(C1 x, C2 y)
        if((isComplex!C1 || isFloatingPoint!C1) && (isComplex!C2 || isFloatingPoint!C2))
{
    static if(isFloatingPoint!C1)
        return .approxEqualC(stdComplex(x, 0), y);
    else
    {
        static if(isFloatingPoint!C2)
            return .approxEqualC(x, stdComplex(y, 0));
        else
        {
            if(x.sqAbs < 1e4) {
                auto e = x - y;
                return std.math.isClose(e.re, 0, 0, 1e-5) && std.math.isClose(e.im, 0, 0, 1e-5);
            } else {
                return std.math.isClose(x.re, y.re) && std.math.isClose(x.im, y.im);
            }
        }
    }
}



typeof(C1.init * C2.init) poly(C1, C2)(C1 x, in C2[] A)
{
    static if(isFloatingPoint!C1 && isFloatingPoint!C2)
        return std.math.poly(x, A);
    else
    {
        // copy from std.math.poly
        immutable size_t N = A.length;
        typeof(return) r = A[N - 1];
        foreach (i; 1 .. N)
        {
            r *= x;
            r += A[N - 1 - i];
        }

        return r;
    }
}
