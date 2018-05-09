module dffdd.window;

import std.math;
import std.traits;
import dffdd.math.math;



/**
点数Nのカイザー窓を作成します
*/
F[] kaiser(F)(size_t N, F beta)
if(isFloatingPoint!F)
{
    F hn = N / F(2);

    immutable F bI0 = besselI0(beta);

    F[] coefs = new F[N];
    foreach(i; 0 .. N){
        immutable F r = F(2)*i / (N-1) - 1;
        immutable F num = sqrt(1-r^^2) * beta;
        // immutable F num = beta * sqrt(1 - ((i - hn)/hn)^^2);
        coefs[i] = besselI0(num) / bI0;
    }

    return coefs;
}

///
unittest
{
    auto window = kaiser!real(12, 14);
    real[] numpyResults = [
        7.72686684e-06, 3.46009194e-03, 4.65200189e-02, 2.29737120e-01,
        5.99885316e-01, 9.45674898e-01, 9.45674898e-01, 5.99885316e-01,
        2.29737120e-01, 4.65200189e-02, 3.46009194e-03, 7.72686684e-06
    ];

    foreach(i; 0 .. window.length)
        assert(approxEqual(window[i], numpyResults[i]));

    window = kaiser!real(13, 14);
    numpyResults = [
        7.72686684e-06, 2.58843902e-03, 3.28855260e-02, 1.64932188e-01,
        4.62716498e-01, 8.28089414e-01, 1.00000000e+00, 8.28089414e-01,
        4.62716498e-01, 1.64932188e-01, 3.28855260e-02, 2.58843902e-03,
        7.72686684e-06];

    foreach(i; 0 .. window.length)
        assert(approxEqual(window[i], numpyResults[i]));
}
