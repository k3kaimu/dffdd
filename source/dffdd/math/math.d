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
        while(abs(r) > F.epsilon){
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
