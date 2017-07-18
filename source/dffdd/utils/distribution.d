module dffdd.utils.ditribution;

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
