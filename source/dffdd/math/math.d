module dffdd.math.math;

import std.math;
import std.traits;


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