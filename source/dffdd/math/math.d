module dffdd.math.math;

import std.math;
import std.traits;

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
