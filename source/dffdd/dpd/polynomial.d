module dffdd.dpd.polynomial;

import std.algorithm;
import std.complex;
import std.math;

import dffdd.utils.linalg;
import dffdd.utils.unit;



enum DPDMode : string
{
    Linearity = "ModeL",
    EfficiencyAndLinearity = "ModeEL"
}



struct PolynomialPredistorter(C)
{
    private {
        alias R = typeof(C.init.re);
    }

    /**
    非線形次数と目標の利得を指定して構築
    */
    this(size_t Porder)
    {
        // _targetGain = targetGain;

        _coefs ~= C(1);
        foreach(i; 1 .. (Porder+1)/2)
            _coefs ~= C(0);
    }


    /**
    PAに入力した信号と，その時の出力信号，及び目標となる出力電力からDPDの係数を推定する．
    目標電力がVoltage(0)に等しいときは，与えられた出力信号から平均電力を算出して，その値を目標にしてDPDの係数を推定する．
    */
    void estimate(in C[] transmitted, in C[] received, Voltage targetPower = Voltage(0))
    in(transmitted.length == received.length)
    do {
        immutable inputPower = transmitted.map!"a.re^^2 + a.im^^2".sum() / received.length;
        immutable outputPower = received.map!"a.re^^2 + a.im^^2".sum() / received.length;

        if(targetPower == Voltage(0))
            targetPower = Voltage.fromWatt(outputPower);

        immutable normCoefIN = sqrt(inputPower);
        immutable normCoefOT = sqrt(targetPower.watt);
        // immutable normCoef = 1;

        auto lsEst = LeastSquareEstimator!(C)(_coefs.length);
        C[] xvec = new C[_coefs.length];
        foreach(i; 0 .. received.length) {
            auto x = received[i] / normCoefOT;
            auto y = transmitted[i] / normCoefIN;
            foreach(p, ref e; xvec)
                e = x * sqAbs(x) ^^ p;
            
            lsEst.append(xvec, y);
        }

        _coefs[] = lsEst.estimate()[];
        _trained = true;
    }


    void apply(in C[] input, ref C[] output)
    {
        if(output.length != input.length)
            output.length = input.length;

        // auto g0 = _targetGain.asV;

        foreach(i, e; input) {
            auto x = e;
            output[i] = C(0);
            foreach(p, c; _coefs)
                output[i] += c * x * x.sqAbs() ^^ p;
        }
    }


    inout(C)[] coefs() inout { return _coefs; }


    bool isConverged()
    {
        return _trained;
    }


  private:
    // Gain _targetGain;
    C[] _coefs;
    bool _trained;
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
        ys ~= C(cbrt(xs[$-1].re));
    }

    auto dist = PolynomialPredistorter!C(5);
    dist.estimate(xs, ys);
    assert(isClose(dist._coefs[0].re, 0, 0, 1e-5));
    assert(isClose(dist._coefs[0].im, 0, 0, 1e-5));
    // assert(isClose(dist._coefs[1].re, 0.843434));
    assert(isClose(dist._coefs[1].im, 0, 0, 1e-5));
    assert(isClose(dist._coefs[2].re, 0, 0, 1e-5));
    assert(isClose(dist._coefs[2].im, 0, 0, 1e-5));

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


    auto gain = (ys[$-1] / ws[$-1]).re;

    foreach(i; 0 .. 100) {
        assert(isClose(ys[i].re, ws[i].re * gain));
        assert(isClose(ys[i].im, ws[i].im * gain));
    }
    
}
