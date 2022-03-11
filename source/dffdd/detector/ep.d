module dffdd.detector.ep;

import std.algorithm : stdmap = map, min, max;
import std.math : sqrt, isClose, SQRT2, SQRT1_2;
import std.numeric : stddot = dotProduct;
import std.traits : isFloatingPoint, isArray;
import std.typecons : Tuple, tuple;

import dffdd.detector.primitives;
import dffdd.mod.primitives : softDecision, softDecision2PAM_BPSK, softDecision2PAM_QPSK, softDecision4PAM_16QAM, softDecision8PAM_64QAM;

// import dffdd.damp.damp;
import dffdd.math.complex;
import dffdd.math.matrix;
import dffdd.math.vector;
import dffdd.math.linalg;
import dffdd.math.exprtemplate : hasMemoryView;
import dffdd.math.matrixspecial : identity;

import dffdd.mod.qpsk;
import dffdd.mod.qam;

import mir.ndslice : slice, sliced, Contiguous, SliceKind;





class EPDectector(C, Mod) : IDetector!(C, C)
if(is(Mod == QPSK!C) || is(Mod == QAM!C) && isComplex!C && (is(typeof(C.init.re) == float) || is(typeof(C.init.re) == double)) )
{
    alias F = typeof(C.init.re);

    alias InputElementType = C;
    alias OutputElementType = C;


    this(Mat)(Mod mod, in Mat chMat, double N0, size_t maxIter = 20)
    if(isMatrixLike!Mat)
    {
        _mod = mod;
        _chMat = matrix!F(chMat.length!0 * 2, chMat.length!1 * 2);
        _rowScales = vector!F(chMat.length!0 * 2, 1);
        _chMat[] = chMat.toRealRIIR;
        _chMatU = matrix!F(chMat.length!0 * 2, chMat.length!0 * 2);
        _chMatV = matrix!F(chMat.length!1 * 2, chMat.length!1 * 2);
        _chMatSigma = vector!F(min(chMat.length!0, chMat.length!1) * 2);

        // // 大システム極限での正規化条件（lim |a_n|^2 = 1）
        // foreach(i; 0 .. chMat.length!0 * 2) {
        //     auto rvec = _chMat.rowVec(i);
        //     auto p = 1/sqrt(norm2(rvec));
        //     _rowScales[i] = p;
        //     rvec.sliced()[] *= p;
        // }

        import kaleidic.lubeck : svd;
        auto svdResult = svd(_chMat.sliced);
        _chMatU[] = svdResult.u.lightScope.matrixed;
        _chMatV[] = svdResult.vt.lightScope.matrixed;
        _chMatSigma[] = svdResult.sigma.lightScope.vectored;

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
            foreach(i; 0 .. inpvecs.length)
                inpvecs[i] *= _rowScales[i];

            // import std.stdio;
            // writeln(inpvecs.sliced);
            Vector!(F, Contiguous) res;

            switch(_mod.symInputLength) {
                case 2: // QPSK
                    res = EP!(softDecision2PAM_QPSK!F, F)
                        (inpvecs, _chMat, _chMatU, _chMatSigma, _chMatV, 0.5, _N0, _maxIter);
                    break;
                case 4: // 16QAM
                    res = EP!(softDecision4PAM_16QAM!F, F)
                        (inpvecs, _chMat, _chMatU, _chMatSigma, _chMatV, 0.5, _N0, _maxIter);
                    break;
                case 6: // 64QAM
                    res = EP!(softDecision8PAM_64QAM!F, F)
                        (inpvecs, _chMat, _chMatU, _chMatSigma, _chMatV, 0.5, _N0, _maxIter);
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
    Vector!(F, Contiguous) _rowScales;
    Matrix!(F, Contiguous) _chMatU, _chMatV;
    Vector!(F, Contiguous) _chMatSigma;
    double _N0;
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

    auto detector = new EPDectector!(C, QPSK!C)(QPSK!C(), chMat, 0.01, 20);
    C[] dst;
    dst = detector.detect(recv, dst);
    // writeln(dst);
    assert(dst[0].re.isClose(recv[0].re));
    assert(dst[0].im.isClose(recv[0].im));
    assert(dst[1].re.isClose(recv[1].re));
    assert(dst[1].im.isClose(recv[1].im));
}


Vector!(F, Contiguous)
    EP(alias prox, F, VecY, MatA)
        (in VecY y_, in MatA A_, in F theta0, in F sigma2, size_t nIteration)
if(isFloatingPoint!F && isVectorLike!VecY && hasMemoryView!VecY && isMatrixLike!MatA && hasMemoryView!MatA && is(VecY.ElementType == F) && is(MatA.ElementType == F))
{
    import kaleidic.lubeck : svd;
    auto svdResult = svd(A_.lightConst.sliced);
    auto chMatU = svdResult.u.lightScope.matrixed;
    auto chMatV = svdResult.v.lightScope.matrixed;
    auto chMatSigma = svdResult.sigma.lightScope.vectored;

    return EP!(prox, F)(y_, A_, chMatU, chMatSigma, chMatV, theta0, sigma2, nIteration);
}


Vector!(F, Contiguous)
    EP(alias prox, F, VecY, MatA, MatAU, VecAS, MatAV)
        (in VecY y_, in MatA A_, in MatAU AU_, in VecAS AS_, in MatAV AV_, in F theta0, in F sigma2, size_t nIteration)
if(isFloatingPoint!F && isVectorLike!VecY && hasMemoryView!VecY && isMatrixLike!MatA && hasMemoryView!MatA && is(VecY.ElementType == F) && is(MatA.ElementType == F))
in(A_.length!0 <= A_.length!1)
{
    import mir.ndslice : map;
    import mir.math.stat : mean;

    auto y = y_.lightConst;
    auto A = A_.lightConst;
    auto U = AU_.lightConst;
    auto V = AV_.lightConst;

    // Mが観測数，Nが未知変数の数
    immutable size_t M = A.length!0,
                     N = A.length!1;

    auto x_AB = vector!F(N, 0);
    auto x_BA = vector!F(N, 0);
    auto x_B = vector!F(N, 0);

    // auto invXi = matrix!F(M, M, 0);
    auto AHinvXi = matrix!F(N, M, 0);

    auto tmpVecM = vector!F(M, 0);
    auto tmpVecN = vector!F(N, 0);
    // auto matAA = A * A.H;

    // auto matAA = matrix!F(M, M, 0);
    // matAA[] = A * A.H;

    F v_AB = 0,
      v_BA = theta0,
      v_B = 0,
      gamma_vBA = 0;

    alias ProxResult = typeof(prox(F(0), theta0));
    auto pxs = slice!ProxResult(N);
    import std.stdio;

    // writeln("START");
    foreach(iter; 1 .. nIteration)
    {
        // SVDの結果を使ってXiの逆行列とA.Hの積とtrace(invXi * matAA)を計算
        F traceInvXiAA = 0;
        AHinvXi[] = F(0);
        foreach(i; 0 .. M) {
            immutable eig = 1/(sigma2 + v_BA * AS_[i]^^2);
            AHinvXi[i, i] = eig * AS_[i];
            traceInvXiAA += eig * AS_[i]^^2;
        }

        tmpVecM[] = y - A * x_BA;
        tmpVecM[] = U.T * tmpVecM;
        tmpVecN[] = AHinvXi * tmpVecM;
        tmpVecN[] = V.T * tmpVecN;

        gamma_vBA = N / traceInvXiAA;
        x_AB[] = x_BA + gamma_vBA * tmpVecN;
        v_AB = gamma_vBA - v_BA;
        
        pxs[] = x_AB.sliced.map!(e => prox(e, v_AB));
        x_B[] = pxs.map!"a.value".vectored;
        v_B = pxs.map!"a.var".mean;

        v_BA = (v_B * v_AB)/(v_AB - v_B);
        x_BA[] = (x_B * v_AB - x_AB * v_B)/(v_AB - v_B);
    }
    // writeln("STOP");

    return x_B;
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

    auto xhat2 = EP!((e, v) => softDecision!(arrP, arrR)(e, v))(y, A, theta0, 0.1 , 20);
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

    auto xhat3 = EP!((e, v) => softDecision!(arrP, arrR)(e, v), F)(y, A, theta0, 0.1 , 200);
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
