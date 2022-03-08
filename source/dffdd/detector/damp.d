module dffdd.detector.damp;

import std.algorithm : stdmap = map;
import std.math : sqrt, isClose, SQRT2, SQRT1_2;
import std.numeric : stddot = dotProduct;
import std.traits : isFloatingPoint, isArray;
import std.typecons : Tuple, tuple;

import dffdd.detector.primitives;
import dffdd.mod.primitives : softDecision;

// import dffdd.damp.damp;
import dffdd.math.complex;
import dffdd.math.matrix;
import dffdd.math.vector;
import dffdd.math.linalg;
import dffdd.math.exprtemplate : hasMemoryView;

import dffdd.mod.qpsk;
import dffdd.mod.qam;

import mir.ndslice : slice, sliced, Contiguous, SliceKind;


class DAMPDectector(C, Mod) : IDetector!(C, C)
if(is(Mod == QPSK!C) || is(Mod == QAM!C) && isComplex!C && (is(typeof(C.init.re) == float) || is(typeof(C.init.re) == double)) )
{
    alias F = typeof(C.init.re);

    alias InputElementType = C;
    alias OutputElementType = C;


    this(Mat)(Mod mod, in Mat chMat, size_t maxIter = 20, F damp1_1 = 1, F damp1_2 = 1, F damp2_1 = 1, F damp2_2 = 1)
    if(isMatrixLike!Mat)
    {
        _mod = mod;
        _chMat = matrix!F(chMat.length!0 * 2, chMat.length!1 * 2);
        _rowScales = vector!F(chMat.length!0 * 2);
        _chMat[] = chMat.toRealRIIR;

        // 大システム極限での正規化条件（lim |a_n|^2 = 1）
        foreach(i; 0 .. chMat.length!0 * 2) {
            auto rvec = _chMat.rowVec(i);
            auto p = 1/sqrt(norm2(rvec));
            _rowScales[i] = p;
            rvec.sliced()[] *= p;
        }

        _delta = _chMat.length!0 * 1.0 / _chMat.length!1;
        _maxIter = maxIter;
        _damp1_1 = damp1_1;
        _damp1_2 = damp1_2;
        _damp2_1 = damp2_1;
        _damp2_2 = damp2_2;
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
            foreach(i; 0 .. inpvecs.length)
                inpvecs[i] *= _rowScales[i];

            // import std.stdio;
            // writeln(inpvecs.sliced);
            Vector!(F, Contiguous) res;

            switch(_mod.symInputLength) {
                case 2: // QPSK
                    res = BODAMP!(softThrBayesOpt2PAM_QPSK!F, F)
                        (inpvecs, _chMat, _delta, 0.5, _maxIter, _damp1_1, _damp1_2, _damp2_1, _damp2_2);
                    break;
                case 4: // 16QAM
                    res = BODAMP!(softThrBayesOpt4PAM_16QAM!F, F)
                        (inpvecs, _chMat, _delta, 0.5, _maxIter, _damp1_1, _damp1_2, _damp2_1, _damp2_2);
                    break;
                case 6: // 64QAM
                    res = BODAMP!(softThrBayesOpt8PAM_64QAM!F, F)
                        (inpvecs, _chMat, _delta, 0.5, _maxIter, _damp1_1, _damp1_2, _damp2_1, _damp2_2);
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


//   private:
    Mod _mod; 
    Matrix!(F, Contiguous) _chMat;
    Vector!(F, Contiguous) _rowScales;
    double _delta;
    double _damp1_1;
    double _damp1_2;
    double _damp2_1;
    double _damp2_2;
    size_t _maxIter;
}

unittest
{
    import std.stdio;
    import mir.ndslice : sliced;
    import dffdd.math.exprtemplate;

    alias C = MirComplex!float;
    auto chMat = [C(1, 0), C(0, 0), C(0, 0), C(1, 0)].sliced(2, 2).matrixed;
    auto recv = [C(1/SQRT2, 1/SQRT2), C(-1/SQRT2, 1/SQRT2)];

    auto detector = new DAMPDectector!(C, QPSK!C)(QPSK!C(), chMat, 20);
    C[] dst;
    dst = detector.detect(recv, dst);
    assert(dst[0].re.isClose(recv[0].re));
    assert(dst[0].im.isClose(recv[0].im));
    assert(dst[1].re.isClose(recv[1].re));
    assert(dst[1].im.isClose(recv[1].im));
}


alias softThrBayesOpt2PAM_BPSK(F) = softDecision!([0.5, 0.5], [-1, 1], F);
alias softThrBayesOpt2PAM_QPSK(F) = softDecision!([0.5, 0.5], [-SQRT1_2, SQRT1_2], F);
alias softThrBayesOpt4PAM_16QAM(F) = softDecision!([0.25, 0.25, 0.25, 0.25], [-3/sqrt(10.0L), -1/sqrt(10.0L), 1/sqrt(10.0L), 3/sqrt(10.0L)], F);
alias softThrBayesOpt8PAM_64QAM(F) = softDecision!([0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125], [-7/sqrt(42.0L), -5/sqrt(42.0L), -3/sqrt(42.0L), -1/sqrt(42.0L), 1/sqrt(42.0L), 3/sqrt(42.0L), 5/sqrt(42.0L), 7/sqrt(42.0L)], F);

// Tuple!(F[], "x_hat", F[], "arrMSE")

Vector!(F, Contiguous)
    BODAMP(alias prox, F, VecY, MatA)
        (in VecY y_, in MatA A_, in F delta, in F theta0, size_t nIteration, F damp1_1, F damp1_2, F damp2_1, F damp2_2)
if(isFloatingPoint!F && isVectorLike!VecY && hasMemoryView!VecY && isMatrixLike!MatA && hasMemoryView!MatA && is(VecY.ElementType == F) && is(MatA.ElementType == F))
{
    auto y = y_.lightConst;
    auto A = A_.lightConst;

    // Mが観測数，Nが未知変数の数
    immutable size_t M = A.length!0,
                     N = A.length!1;

    auto x = vector!F(N, 0);
    auto z = vector!F(M, 0);
    F theta = theta0;
    auto w = vector!F(N, 0);

    import mir.ndslice : map, slice;
    import mir.math.stat : mean;
    import dffdd.math.linalg : norm2;

    alias ProxResult = typeof(prox(F(0), theta0^^2/delta));
    auto wprox = slice!ProxResult(N);
    wprox[] = w.sliced.map!(e => prox(e, theta^^2/delta));

    foreach(iter; 1 .. nIteration) {
        F damp2 = damp2_2 * iter / nIteration + damp2_1 * (1 - 1.0 * iter / nIteration);
        F damp1 = damp1_2 * iter / nIteration + damp1_1 * (1 - 1.0 * iter / nIteration);

        z[] = (1-damp1) * z + damp1 * (y - A * x + (1/delta) * z * wprox.map!"a.drv".mean);
        theta = sqrt(norm2(z) / N);
        w[] = x + A.T * z;
        wprox[] = w.sliced.map!(e => prox(e, theta^^2/delta));
        x[] = (1-damp2) * x + damp2 * wprox.map!"a.value".vectored;
    }

    return x;
}


unittest
{
    import std.stdio;


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
    enum F[] arrP = [0.5, 0, 0.5];
    enum F[] arrR = [-1, 0, 1];
    auto theta0 = sqrt(stddot(arrP, arrR.stdmap!"a^^2"));

    int nIteration = 2;
    import std.stdio;

    auto xhat1 = BODAMP!((e, v) => softDecision!(arrP, arrR)(e, v))(y, A, delta, theta0, 3, 1, 1, 1, 1);
    assert(xhat1[0].isClose(-0.9729, 1e-3));
    assert(xhat1[1].isClose(+0.9744, 1e-3));
    assert(xhat1[2].isClose(-0.9220, 1e-3));
    assert(xhat1[3].isClose(+0.9716, 1e-3));
    assert(xhat1[4].isClose(+0.8439, 1e-3));
    assert(xhat1[5].isClose(-0.9594, 1e-3));
    assert(xhat1[6].isClose(+0.9440, 1e-3));
    assert(xhat1[7].isClose(-0.8160, 1e-3));
    assert(xhat1[8].isClose(-0.9436, 1e-3));
    assert(xhat1[9].isClose(+0.8961, 1e-3));

    auto xhat2 = BODAMP!((e, v) => softDecision!(arrP, arrR)(e, v))(y, A, delta, theta0, 20, 1, 1, 1, 1);
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

    auto xhat3 = BODAMP!((e, v) => softDecision!(arrP, arrR)(e, v), F)(y, A, delta, theta0, 200, 1, 1, 1, 1);
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


// unittest
// {
//     import dffdd.math.exprtemplate;
//     import std.stdio;
//     import mir.ndslice;

//     alias F = float;

//     auto chMat = matrix!F(4, 4);
//     chMat[] = ConstEye!(F, 1)(4);
//     writeln(chMat);

//     enum F[] arrR = [-SQRT1_2, SQRT1_2];

//     auto vec = (cast(F[])[SQRT1_2, -SQRT1_2, SQRT1_2, SQRT1_2]).sliced.vectored;
//     // pragma(msg, typeof(vec));
//     // auto res = BODAMP!(softThrBayesOpt!([0.5, 0.5], arrR, F), F)(vec, chMat, 1, 0.5, 50);
//     auto res = BODAMP!(softThrBayesOpt2PAM_QPSK!F, F)(vec, chMat, 1, 0.5, 50);

//     writeln(res);
//     writeln(makeViewOrNewSlice(vec - chMat * res, matvecAllocator));
// }
