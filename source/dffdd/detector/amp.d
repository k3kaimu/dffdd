module dffdd.detector.amp;

import std.math;
import std.traits;

import dffdd.math.complex;
import dffdd.math.matrix;
import dffdd.math.vector;
import dffdd.math.exprtemplate;
import dffdd.math.linalg : toRealRI, fromRealRI, toRealRIIR;

import dffdd.detector.primitives;
import dffdd.mod.primitives;
import dffdd.mod.qpsk;
import dffdd.mod.qam;

import mir.ndslice : Contiguous, sliced;


final class AMPDetector(C, Mod) : IDetector!(C, C)
if((is(Mod == QPSK!C) || is(Mod == QAM!C)) && isComplex!C && (is(typeof(C.init.re) == float) || is(typeof(C.init.re) == double)))
{
    alias F = typeof(C.init.re);
    alias InputElementType = C;
    alias OutputElementType = C;


    this(Mat)(Mod mod, in Mat chMat, double N0, size_t maxIter = 20)
    if(isMatrixLike!Mat)
    {
        _mod = mod;
        _chMat = matrix!F(chMat.length!0 * 2, chMat.length!1 * 2);
        _chMat[] = chMat.toRealRIIR;

        _delta = _chMat.length!0 * 1.0 / _chMat.length!1;
        _maxIter = maxIter;
        _N0 = N0;
    }


    size_t inputLength() const { return _chMat.length!0 / 2; }
    size_t outputLength() const { return _chMat.length!1 / 2; }


    C[] detect(in C[] received, return ref C[] detected) const
    in(received.length % this.inputLength == 0)
    {
        immutable M = this.inputLength;
        immutable N = this.outputLength;

        if(detected.length != received.length / M * N)
            detected.length = received.length / M * N;

        // auto reals = received.toRealRI;
        auto inpvecs = vector!F(M * 2);
        foreach(n; 0 .. received.length / M) {
            inpvecs[] = received[M * n  .. M * (n+1)].sliced.vectored.toRealRI;
            // foreach(i; 0 .. inpvecs.length)
            //     inpvecs[i] *= _rowScales[i];

            Vector!(F, Contiguous) res;
            switch(_mod.symInputLength) {
                case 2: // QPSK
                    res = AMP!(softDecision2PAM_QPSK!F, F)
                        (inpvecs, _chMat, _delta, 0.5, _N0, _maxIter);
                    break;
                case 4: // 16QAM
                    res = AMP!(softDecision4PAM_16QAM!F, F)
                        (inpvecs, _chMat, _delta, 0.5, _N0, _maxIter);
                    break;
                case 6: // 64QAM
                    res = AMP!(softDecision8PAM_64QAM!F, F)
                        (inpvecs, _chMat, _delta, 0.5, _N0, _maxIter);
                    break;
                default: {
                    import std.conv;
                    assert(0, "Unsupported modulation scheme: " ~ _mod.to!string);
                }
            }
            // writeln(res);
            detected[N * n .. N * (n+1)].sliced.vectored[] = res.fromRealRI;
        }

        return detected;
    }


  private:
    Mod _mod;
    Matrix!(F, Contiguous) _chMat;
    // Vector!(F, Contiguous) _rowScales;
    double _delta;
    size_t _maxIter;
    double _N0;
}


unittest
{
    import std.stdio;
    import mir.ndslice : sliced;
    import dffdd.math.exprtemplate;

    alias C = MirComplex!float;
    auto chMat = [C(1, 0), C(0, 0), C(0, 0), C(1, 0)].sliced(2, 2).matrixed;
    auto recv = [C(1/SQRT2, 1/SQRT2), C(-1/SQRT2, 1/SQRT2)];

    auto detector = new AMPDetector!(C, QPSK!C)(QPSK!C(), chMat, 0.1, 20);
    C[] dst;
    dst = detector.detect(recv, dst);
    assert(dst[0].re.isClose(recv[0].re));
    assert(dst[0].im.isClose(recv[0].im));
    assert(dst[1].re.isClose(recv[1].re));
    assert(dst[1].im.isClose(recv[1].im));
}


Vector!(F, Contiguous)
    AMP(alias prox, F, VecY, MatA)
        (in VecY y_, in MatA A_, in F delta, in F theta0, in F sigma2, size_t nIteration)
if(isFloatingPoint!F && isVectorLike!VecY && hasMemoryView!VecY && isMatrixLike!MatA && hasMemoryView!MatA && is(VecY.ElementType == F) && is(MatA.ElementType == F))
{
    auto y = y_.lightConst;
    auto A = A_.lightConst;

    // Mが観測数，Nが未知変数の数
    immutable size_t M = A.length!0,
                     N = A.length!1;

    auto x_AB = vector!F(N, 0);
    auto x_B = vector!F(N, 0);
    auto u = vector!F(M, 0);

    F v_B = theta0;
    F v_AB = sigma2 + (1/delta) * v_B;

    import mir.ndslice : map, slice;
    import mir.math.stat : mean;

    alias ProxResult = typeof(prox(F(0), theta0));
    auto xprox = slice!ProxResult(N);
    // xprox[] = w.sliced.map!(e => prox(e, theta^^2/delta));
    import std.stdio;
    foreach(iter; 1 .. nIteration) {
        u[] = y - A * x_B + (v_B/delta/v_AB) * u;
        x_AB[] = x_B + A.T * u;
        v_AB = sigma2 + (1/delta) * v_B;
        xprox[] = x_AB.sliced.map!(e => prox(e, v_AB));
        x_B[] = xprox.map!"a.value".vectored;
        v_B = xprox.map!"a.var".mean;
    }

    return x_B;
}


unittest
{
    import std.stdio;
    import std.algorithm : stdmap = map;
    import std.numeric : stddot = dotProduct;
    import mir.ndslice;
    alias F = double;

    auto A = (cast(F[])[0.154508497187474,  -0.404508497187474,   0.500000000000000, -0.404508497187474,  0.154508497187474,    0.154508497187474,  -0.404508497187474,   0.500000000000000, -0.404508497187473,  0.154508497187473,
    0.226995249869773,  -0.493844170297569,   0.353553390593274, 0.0782172325201156,  -0.445503262094184,   0.445503262094184, -0.0782172325201147, -0.353553390593274,  0.493844170297569,    -0.226995249869773,
    0.445503262094184,  0.0782172325201155,   -0.353553390593274,    -0.493844170297569, -0.226995249869773,  0.226995249869773,    0.493844170297569,  0.353553390593274,    -0.0782172325201148,    -0.445503262094184,
    0.475528258147577,  0.293892626146237,    3.06161699786838e-17,   -0.293892626146237,    -0.475528258147577, -0.475528258147577,  -0.293892626146237,   -9.18485099360515e-17, 0.293892626146236,   0.475528258147577,
    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,
    0.293892626146237,  -0.475528258147577,   -9.18485099360515e-17, 0.475528258147577,   -0.293892626146236,    -0.293892626146237, 0.475528258147577,   2.75545529808155e-16,  -0.475528258147577,   0.293892626146237,
    0.493844170297569,  0.445503262094184,    0.353553390593274,  0.226995249869773,    0.0782172325201155, -0.0782172325201154, -0.226995249869773,  -0.353553390593274,   -0.445503262094184,    -0.493844170297569,
    0.353553390593274,  -0.353553390593274,   -0.353553390593274,    0.353553390593274,  0.353553390593274,    -0.353553390593273, -0.353553390593274,  0.353553390593273,    0.353553390593274,  -0.353553390593273,
    ]).sliced(8, 10).matrixed;

    auto y = (cast(F[])[
        -1.81199936706226,
        -2.35139170510204,
        -1.10222286082885,
        -0.555771890670480,
        0.0283159355873546,
        0.911307827117349,
        0.0769074757814163,
        -0.686546300398577,]).sliced.vectored;

    auto delta = 0.8;
    enum F[] arrP = [0.5, 0.5];
    enum F[] arrR = [-1, 1];
    auto theta0 = sqrt(stddot(arrP, arrR.stdmap!"a^^2"));

    int nIteration = 2;
    import std.stdio;

    auto xhat2 = AMP!((e, v) => softDecision!(arrP, arrR)(e, v))(y, A, delta, theta0, 0.1, 20);
    // writeln(xhat2);
    assert(xhat2[0].isClose(-1, 1e-3));
    assert(xhat2[1].isClose(+1, 1e-3));
    assert(xhat2[2].isClose(-1, 1e-3));
    assert(xhat2[3].isClose(+1, 1e-3));
    assert(xhat2[4].isClose(+1, 1e-3));
    assert(xhat2[5].isClose(-1, 1e-3));
    assert(xhat2[6].isClose(+1, 1e-3));
    assert(xhat2[7].isClose(-1, 1e-3));
    assert(xhat2[8].isClose(-1, 1e-3));
    assert(xhat2[9].isClose(+1, 1e-3));

    auto xhat3 = AMP!((e, v) => softDecision!(arrP, arrR)(e, v), F)(y, A, delta, theta0, 0.1, 200);
    // writeln(xhat3);
    assert(xhat3[0].isClose(-1, 1e-3));
    assert(xhat3[1].isClose(+1, 1e-3));
    assert(xhat3[2].isClose(-1, 1e-3));
    assert(xhat3[3].isClose(+1, 1e-3));
    assert(xhat3[4].isClose(+1, 1e-3));
    assert(xhat3[5].isClose(-1, 1e-3));
    assert(xhat3[6].isClose(+1, 1e-3));
    assert(xhat3[7].isClose(-1, 1e-3));
    assert(xhat3[8].isClose(-1, 1e-3));
    assert(xhat3[9].isClose(+1, 1e-3));
}

