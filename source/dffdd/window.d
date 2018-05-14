module dffdd.window;

import std.algorithm;
import std.math;
import std.traits;
import dffdd.math.math;


/**
阻止帯域での減衰量からカイザー窓のパラメータbetaを算出します
See also: https://jp.mathworks.com/help/signal/ref/kaiser.html
*/
real kaiserBetaFromStopbandAttdB(real attdB) pure nothrow @safe @nogc
{
    real beta;
    if(attdB >= 50)
        beta = 0.1102 * (attdB - 8.7);
    else if(attdB > 21)
        beta = 0.5842 * (attdB - 21)^^0.4 + 0.07886 * (attdB - 21);
    else
        beta = 0;

    return beta;
}

unittest
{
    assert(approxEqual(kaiserBetaFromStopbandAttdB(0), 0));
    assert(approxEqual(kaiserBetaFromStopbandAttdB(22), 0.5842 + 0.07886));
    assert(approxEqual(kaiserBetaFromStopbandAttdB(50), 4.55126));
}


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


/**
点数Nのカイザー・ベッセル派生窓を作成します
*/
F[] kaiserBesselDerived(F)(size_t N, F beta)
in{
    assert(N % 2 == 0);
}
do {
    immutable size_t M = N / 2;

    auto kwin = kaiser!F(M+1, beta);
    auto dst = new F[N];

    immutable F sumk = kwin.sum();

    F partsum = 0;
    foreach(i; 0 .. N/2){
        partsum += kwin[i];
        dst[i] = sqrt(partsum / sumk);
        dst[$ - 1 - i] = dst[i];
    }


    return dst;
}
