module dffdd.utils.distribution;

import std.complex;
import std.math;
import std.random;
import dffdd.math.math;

F normalDist(F = real, Rnd)(ref Rnd rnd)
{
    F x = uniform01(rnd),
      y = uniform01(rnd);

    return sqrt(-2 * log(x)) * cos(2 * PI * y);
}


F normalDist(F, Rnd)(F mu, F sigma, ref Rnd rnd)
{
    auto x = normalDist(rnd);
    return x * sigma + mu;
}


Complex!F complexGaussian01(F = real, Rnd)(ref Rnd rnd)
{
    F x = uniform01(rnd),
      y = uniform01(rnd);

    typeof(return) dst = fast_sqrt!F(-1 * fast_log!F(x)) * fast_expi!(Complex!F)(2 * PI * y);

    return dst;
}


Complex!F complexGaussian(F = real, Rnd)(Complex!F mu, F sigma, ref Rnd rnd)
{
    auto x = complexGaussian01!F(rnd);
    return x * sigma + mu;
}
