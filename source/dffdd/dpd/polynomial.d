module dffdd.dpd.polynomial;

import std.algorithm;
import std.complex;
import std.math;

import dffdd.utils.linalg;
import dffdd.utils.unit;


struct PolynomialPredistorter(C)
{
    /**
    非線形次数と目標の利得を指定して構築
    */
    this(size_t Porder, Gain targetGain)
    {
        _targetGain = targetGain;

        _coefs ~= C(1);
        foreach(i; 1 .. (Porder+1)/2)
            _coefs ~= C(0);
    }


    void estimate(in C[] transmitted, in C[] received)
    in(transmitted.length == received.length)
    do {
        immutable outputPower = received.map!"a.re^^2 + a.im^^2".sum() / received.length;
        immutable normCoef = sqrt(outputPower);

        auto lsEst = LeastSquareEstimator!(C)(_coefs.length);
        C[] xvec = new C[_coefs.length];
        foreach(i; 0 .. received.length) {
            auto x = received[i] / normCoef;
            auto y = transmitted[i];
            foreach(p, ref e; xvec)
                e = x * sqAbs(x) ^^ p;
            
            lsEst.append(xvec, y);
        }

        C[] estparam = lsEst.estimate();
        foreach(p, e; estparam)
            _coefs[p] = estparam[p] / (normCoef ^^ (2*p + 1));
    }


    void apply(in C[] input, ref C[] output)
    {
        if(output.length != input.length)
            output.length = input.length;

        foreach(i, e; input) {
            auto x = e * _targetGain.asV;
            output[i] = C(0);
            foreach(p, c; _coefs)
                output[i] += c * x * x.sqAbs() ^^ p;
        }
    }


  private:
    Gain _targetGain;
    C[] _coefs;
}

unittest
{
    import std.random;
    import std.stdio;

    alias C = Complex!double;

    Random rnd = Random(0);
    C[] xs, ys;
    foreach(i; 0 .. 100) {
        xs ~= C(uniform01(rnd));
        ys ~= C(cbrt(xs[$-1].re)) * 3;
    }

    auto dist = PolynomialPredistorter!C(5, 20.dB);
    dist.estimate(xs, ys);
    assert(approxEqual(dist._coefs[0].re, 0));
    assert(approxEqual(dist._coefs[0].im, 0));
    assert(approxEqual(dist._coefs[1].re, 1/27.0));
    assert(approxEqual(dist._coefs[1].im, 0));
    assert(approxEqual(dist._coefs[2].re, 0));
    assert(approxEqual(dist._coefs[2].im, 0));

    C[] ws;
    xs = [];
    ys = [];
    foreach(i; 0 .. 100) {
        ws ~= C(uniform01(rnd));
    }
    dist.apply(ws, xs);
    foreach(i; 0 .. 100) {
        ys ~= C(cbrt(xs[i].re)) * 3;
    }
    foreach(i; 0 .. 100) {
        assert(approxEqual(ys[i].re, ws[i].re * 10));
        assert(approxEqual(ys[i].im, ws[i].im * 10));
    }
    
}
