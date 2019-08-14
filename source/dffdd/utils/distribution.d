module dffdd.utils.distribution;

import std.math;
import std.random;

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

    return sqrt(-2 * log(x)) * std.complex.expi(2 * PI * y);
}


Complex!F complexGaussian(F = real, Rnd)(Complex!F mu, F sigma, ref Rnd rnd)
{
    auto x = complexGaussian01!F(rnd);
    return x * sigma + mu;
}
