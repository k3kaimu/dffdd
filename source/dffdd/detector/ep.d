module dffdd.detector.ep;

import std.algorithm : stdmap = map, min, max;
import std.math : sqrt, isClose, SQRT2, SQRT1_2;
import std.numeric : stddot = dotProduct;
import std.traits : isFloatingPoint, isArray;
import std.typecons : Tuple, tuple;

import dffdd.detector.primitives;
import dffdd.mod.primitives/* : softDecision, softDecision2PAM_BPSK, softDecision2PAM_QPSK, softDecision4PAM_16QAM, softDecision8PAM_64QAM*/;

// import dffdd.damp.damp;
import dffdd.math.complex;
import dffdd.math.matrix;
import dffdd.math.vector;
import dffdd.math.linalg;
import dffdd.math.exprtemplate : hasMemoryView;
import dffdd.math.matrixspecial : identity, diag;
import dffdd.math.math;

import dffdd.mod.qpsk;
import dffdd.mod.qam;

import mir.ndslice : slice, sliced, Contiguous, SliceKind;



version(LDC)
{
    import ldc.attributes : fastmath;
}
else
{
    enum fastmath = 0;
}



class EPDetector(C, Mod, F = typeof(C.init.re), MatU = Matrix!(F, Contiguous), MatV = Matrix!(F, Contiguous)) : IDetector!(C, C)
if((is(Mod : QPSK!C) || is(Mod : QAM!C)) && isComplex!C && (is(typeof(C.init.re) == float) || is(typeof(C.init.re) == double)) )
{
    alias InputElementType = C;
    alias OutputElementType = C;


  static if(is(MatU == Matrix!(F, Contiguous)) && is(MatV == Matrix!(F, Contiguous)))
  {
    this(Mat)(Mod mod, in Mat chMat, double N0, size_t maxIter = 20)
    if(isMatrixLike!Mat)
    {
        _mod = mod;

        static if(is(F == typeof(C.init.re)))
        {
            auto chm = matrix!F(chMat.length!0 * 2, chMat.length!1 * 2);
            chm[] = chMat.toRealRIIR;

            _rowScales = vector!F(chMat.length!0 * 2, F(1));
            _chMatU = matrix!F(chMat.length!0 * 2, chMat.length!0 * 2);
            _chMatV = matrix!F(chMat.length!1 * 2, chMat.length!1 * 2);
            _chMatSigma = vector!F(min(chMat.length!0, chMat.length!1) * 2);

            import kaleidic.lubeck : svd;
            auto svdResult = svd(chm.sliced);
            _chMatU[] = svdResult.u.lightScope.matrixed;
            _chMatV[] = svdResult.vt.lightScope.matrixed;
            _chMatSigma[] = svdResult.sigma.lightScope.vectored;
        }
        else
        {
            auto chm = matrix!F(chMat.length!0, chMat.length!1);
            chm[] = chMat;
            _rowScales = vector!F(chMat.length!0, F(1));

            _chMatU = matrix!F(chMat.length!0, chMat.length!0);
            _chMatV = matrix!F(chMat.length!1, chMat.length!1);
            _chMatSigma = vector!F(min(chMat.length!0, chMat.length!1));

            import kaleidic.lubeck : svd;
            auto svdResult = svd(chm.sliced);
            _chMatU[] = svdResult.u.lightScope.matrixed;
            _chMatV[] = svdResult.vt.lightScope.matrixed;
            _chMatSigma[] = svdResult.sigma.lightScope.vectored;
        }

        _maxIter = maxIter;
        _N0 = N0;
    }
  }


    this(VecS)(Mod mod, MatU chMatU, VecS chMatSigma, MatV chMatVh, double N0, size_t maxIter = 20)
    {
        _mod = mod;
        _rowScales = vector!F(chMatU.length!0, F(1));
        _chMatU = chMatU;
        _chMatSigma = vector!F(min(chMatU.length!0, chMatVh.length!1));
        _chMatSigma[] = chMatSigma;
        _chMatV = chMatVh;
        _maxIter = maxIter;
        _N0 = N0;
    }


    size_t inputLength() const
    {
        return is(F == typeof(C.init.re)) ? _chMatU.length!0 / 2 : _chMatU.length!0;
    }


    size_t outputLength() const
    {
        return is(F == typeof(C.init.re)) ? _chMatV.length!1 / 2 : _chMatV.length!1;
    }



    C[] detect(in C[] received, return ref C[] detected) /*const*/
    in(received.length % this.inputLength == 0)
    {
        immutable M = this.inputLength;
        immutable N = this.outputLength;

        if(detected.length != received.length / M * N)
            detected.length = received.length / M * N;

        auto inpvecs = vector!F(_chMatU.length!0);
        foreach(n; 0 .. received.length / M) {

            static if(is(F == typeof(C.init.re)))
                inpvecs[] = received[M * n  .. M * (n+1)].sliced.vectored.toRealRI;
            else
                inpvecs[] = received[M * n  .. M * (n+1)].sliced.vectored;

            foreach(i; 0 .. inpvecs.length)
                inpvecs[i] *= _rowScales[i];

            Vector!(F, Contiguous) res;

            switch(_mod.symInputLength) {
                case 2: // QPSK
                    static if(is(F == typeof(C.init.re)))
                    {
                        res = EP!(softDecision2PAM_QPSK!F, F, F)
                            (inpvecs, _chMatU, _chMatSigma, _chMatV, 0.5, _N0, _maxIter);
                    }
                    else
                    {
                        import std.math : SQRT1_2;
                        res = EP!(/*softDecisionQPSK!C*/approx_softDecisionQPSK, C, typeof(C.init.re))
                            (inpvecs, _chMatU, _chMatSigma, _chMatV, 1, _N0, _maxIter);
                    }
                    break;
                case 4: // 16QAM
                    static if(is(F == typeof(C.init.re)))
                    {
                        res = EP!(softDecision4PAM_16QAM!F, F, F)
                            (inpvecs, _chMatU, _chMatSigma, _chMatV, 0.5, _N0, _maxIter);
                    }
                    else
                    {
                        res = EP!(/*softDecision16QAM!C*/approx_softDecision16QAM, C, typeof(C.init.re))
                            (inpvecs, _chMatU, _chMatSigma, _chMatV, 1, _N0, _maxIter);
                    }
                    break;
                case 6: // 64QAM
                    static if(is(F == typeof(C.init.re)))
                    {
                        res = EP!(softDecision8PAM_64QAM!F, F, F)
                            (inpvecs, _chMatU, _chMatSigma, _chMatV, 0.5, _N0, _maxIter);
                    }
                    else
                    {
                        res = EP!(/*softDecision64QAM!C*/approx_softDecision64QAM, C, typeof(C.init.re))
                            (inpvecs, _chMatU, _chMatSigma, _chMatV, 1, _N0, _maxIter);
                    }
                    break;
                default: {
                    import std.conv;
                    assert(0, "Unsupported modulation scheme: " ~ _mod.to!string);
                }
            }
            // writeln(res);

            static if(is(F == typeof(C.init.re)))
                detected[N * n .. N * (n+1)].sliced.vectored[] = res.fromRealRI;
            else
                detected[N * n .. N * (n+1)].sliced.vectored[] = res;
        }

        return detected;
    }


  private:
    Mod _mod; 
    Vector!(F, Contiguous) _rowScales;
    MatU _chMatU;
    MatV _chMatV;
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

    auto detector = new EPDetector!(C, QPSK!C)(QPSK!C(), chMat, 0.01, 20);
    C[] dst;
    dst = detector.detect(recv, dst);
    // writeln(dst);
    assert(dst[0].re.isClose(recv[0].re));
    assert(dst[0].im.isClose(recv[0].im));
    assert(dst[1].re.isClose(recv[1].re));
    assert(dst[1].im.isClose(recv[1].im));
}

unittest
{
    import std.stdio;
    import mir.ndslice : sliced;
    import dffdd.math.exprtemplate;

    alias C = MirComplex!float;
    auto chMat = [C(1, 0), C(0, 0), C(0, 0), C(1, 0)].sliced(2, 2).matrixed;
    auto recv = [C(1/SQRT2, 1/SQRT2), C(-1/SQRT2, 1/SQRT2)];

    auto detector = new EPDetector!(C, QPSK!C, C)(QPSK!C(), chMat, 0.01, 20);
    C[] dst;
    dst = detector.detect(recv, dst);
    // writeln(dst);
    assert(dst[0].re.isClose(recv[0].re));
    assert(dst[0].im.isClose(recv[0].im));
    assert(dst[1].re.isClose(recv[1].re));
    assert(dst[1].im.isClose(recv[1].im));
}


auto makeSVDEPDetector(C, Mod, MatU, VecS, MatV)(Mod mod, MatU matU, VecS vecS, MatV matV, double N0, size_t maxIter = 20)
{
    return new EPDetector!(C, Mod, MatU.ElementType, MatU, MatV)(mod, matU, vecS, matV, N0, maxIter);
}


Vector!(F, Contiguous)
    EP(alias prox, C, F, VecY, MatA)
        (in VecY y_, in MatA A_, in F theta0, in F sigma2, size_t nIteration)
if((isComplex!C || isFloatingPoint!C) && isFloatingPoint!F && isVectorLike!VecY && hasMemoryView!VecY && isMatrixLike!MatA && hasMemoryView!MatA && is(VecY.ElementType == C) && is(MatA.ElementType == C))
{
    import kaleidic.lubeck : svd;
    auto svdResult = svd(A_.lightConst.sliced);
    auto chMatU = svdResult.u.lightScope.matrixed;
    auto chMatV = svdResult.vt.lightScope.matrixed;
    auto chMatSigma = svdResult.sigma.lightScope.vectored;

    return EP!(prox, C, F)(y_, chMatU, chMatSigma, chMatV, theta0, sigma2, nIteration);
}


@fastmath
Vector!(C, Contiguous)
    EP(alias prox, C = VecY.ElementType, F, VecY, MatAU, VecAS, MatAV)
        (VecY y, MatAU U, VecAS AS_, MatAV V, in F theta0, in F sigma2, size_t nIteration)
if((isNarrowComplex!C || isFloatingPoint!C) && isFloatingPoint!F
    && isVectorLike!VecY && hasMemoryView!VecY
    && isMatrixLike!MatAU && isMatrixLike!MatAV && isVectorLike!VecAS
    && is(VecY.ElementType == C) && is(MatAU.ElementType == C) && is(MatAV.ElementType == C)
    && (is(VecAS.ElementType == C) || is(VecAS.ElementType == F)))
in(U.length!0 <= V.length!1)
{
    import std.typecons : scoped;
    import dffdd.math.complex : abs;
    import mir.ndslice : map;
    import mir.math.stat : mean;
    import dffdd.math.exprtemplate : TempMemoryManager;

    // Mが観測数，Nが未知変数の数
    immutable size_t M = U.length!0,
                     N = V.length!0;

    auto tmm = scoped!TempMemoryManager(M * C.sizeof * 8 + N * C.sizeof * 8);

    auto x_AB = tmm.makeVectorFromTMM_SIMD!(C, 4)(N, C(0));
    auto x_BA = tmm.makeVectorFromTMM_SIMD!(C, 4)(N, C(0));
    auto x_B = tmm.makeVectorFromTMM_SIMD!(C, 4)(N, C(0));

    // auto vecAHinvXi = tmm.makeSlice!C(M);
    auto vecAHinvXi = tmm.makeVectorFromTMM_SIMD!(C, 4)(N, C(0));

    auto tmpVecM = tmm.makeVectorFromTMM_SIMD!(C, 4)(M, C(0));
    auto tmpVecN = tmm.makeVectorFromTMM_SIMD!(C, 4)(N, C(0));
    auto tmpVecN2 = tmm.makeVectorFromTMM_SIMD!(C, 4)(N, C(0));

    auto y2 = tmm.makeVectorFromTMM_SIMD!(C, 4)(M, C(0));
    y2.noalias = U.H * y;

    // UTA = U.T * A 
    auto UTA = diag(M, N, AS_.sliced) * V;

    F v_AB = 0,
      v_BA = theta0,
      v_B = 0,
      gamma_vBA = 0;

    enum F EPS = 1e-6;

    // writeln("START");
    foreach(iter; 1 .. nIteration)
    {
        // SVDの結果を使ってXiの逆行列とA.Hの積とtrace(invXi * matAA)を計算
        F traceInvXiAA = EPS;
        vecAHinvXi[] = C(0);
        foreach(i; 0 .. M) {
            immutable ASP = sqAbs(AS_[i]);
            immutable eig = 1/(sigma2 + v_BA * ASP + EPS);
            vecAHinvXi[i] = eig * conj(AS_[i]);
            traceInvXiAA += eig * ASP;
        }

        tmpVecM.noalias = y2 - UTA * x_BA;
        tmpVecN2[] = F(0);
        foreach(i; 0 .. M) tmpVecN2[i] = vecAHinvXi[i] * tmpVecM[i];
        // tmpVecN2.sliced[0 .. M] = tmpVecM.sliced;
        // tmpVecN2.sliced[0 .. M] *= vecAHinvXi;
        tmpVecN.noalias = V.H * tmpVecN2;

        gamma_vBA = N / traceInvXiAA;
        x_AB.noalias = x_BA + gamma_vBA * tmpVecN;
        v_AB = gamma_vBA - v_BA;
        
        // pxs[] = x_AB.sliced.map!(e => prox(e, v_AB));
        // x_B.noalias = pxs.map!"a.value".vectored;
        // v_B = pxs.map!"a.var".mean;
        import core.simd;
        static if(isNarrowComplex!C && is(typeof((float8 re, float8 im, float v){ auto a = prox(re, im, v); })))
        {
            v_B = 0;
            immutable size_t LEN = x_AB.length;
            size_t k = 0;
            while(LEN - k >= 8) {
                float8 re, im;
                static foreach(i; 0 .. 8) {
                    re[i] = x_AB[k + i].re;
                    im[i] = x_AB[k + i].im;
                }

                auto p = prox(re, im, v_AB);
                static foreach(i; 0 .. 8) {
                    x_B[k + i].re = p.re[i]; 
                    x_B[k + i].im = p.im[i];
                    v_B += p.var;
                }

                k += 8;
            }
            {
                float8 re = 0, im = 0;
                foreach(i; 0 .. (LEN - k)) {
                    re[i] = x_AB[k + i].re;
                    im[i] = x_AB[k + i].im;
                }
                auto p = prox(re, im, v_AB);
                foreach(i; 0 .. (LEN - k)) {
                    x_B[k + i].re = p.re[i]; 
                    x_B[k + i].im = p.im[i];
                    v_B += p.var;
                }
            }

            v_B /= x_AB.length;
        }
        else
        {
            v_B = 0;
            foreach(i; 0 .. x_AB.length) {
                auto p = prox(x_AB[i], v_AB);
                x_B[i] = p.value;
                v_B += p.var;
            }
            v_B /= x_AB.length;
        }

        F vden = v_AB - v_B;
        if(fast_abs!F(vden) < EPS)
            vden += EPS;

        v_BA = (v_B * v_AB)/vden;
        x_BA.noalias = (x_B * v_AB - x_AB * v_B)/vden;
    }
    // writeln("STOP");

    // return x_B;

    auto dst = vector!C(N);
    dst[] = x_B;
    return dst;
}


unittest
{
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

    auto xhat2 = EP!((e, v) => softDecision!(arrP, arrR)(e, v), F, F)(y, A, theta0, 0.1, 20);
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

    auto xhat3 = EP!((e, v) => softDecision!(arrP, arrR)(e, v), F, F)(y, A, theta0, 0.1 , 200);
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

unittest
{
    import std.range;
    import std.array;
    import std.math;
    import std.complex;
    import dffdd.math.complex;
    import std.random;
    import dffdd.utils.distribution;
    import dffdd.mod.primitives;

    alias F = double;
    alias C = Complex!double;

    Random rnd;
    rnd.seed(0);

    size_t M = 8;
    size_t N = 10;

    auto W = matrix!C(N, N);
    foreach(i; 0 .. N)
        foreach(j; 0 .. N)
            W[i, j] = std.complex.expi(2*PI/N*i*j) / sqrt(N*1.0);
    
    auto P = matrix!C(N, N, C(0));
    auto idxs = iota(N).array();
    std.random.randomShuffle(idxs, rnd);
    foreach(i; 0 .. N)
        P[i, idxs[i]] = C(1);

    auto V = matrix!C(N, N);
    V[] = P * W;

    auto U = matrix!C(M, M);
    U[] = identity!C(M);

    auto ds = vector!C(M);
    foreach(i; 0 .. M)
        ds[i] = (0.1)^^(i*1.0/M) * std.complex.expi(2*PI*i / M);

    auto D = diag(M, N, ds.sliced);
    auto A = matrix!C(M, N);
    A[] = U * D * V;

    auto x0 = vector!C(N);
    foreach(i; 0 .. N)
        x0[i] = C(choice([-1, 1], rnd), choice([-1, 1], rnd)) * SQRT1_2;

    auto y0 = vector!C(M);
    y0[] = A * x0;

    immutable F sigma2 = 0.01;
    foreach(i; 0 .. M)
        y0[i] += complexGaussian01!F(rnd) * sqrt(sigma2);

    auto xhat1 = EP!(softDecisionQPSK!C, C, F)(y0, U, ds, V, 1, sigma2, 2);

    auto Ar = matrix!F(M * 2, N * 2);
    Ar[] = A.toRealRIIR;

    auto y0r = vector!F(M * 2);
    y0r[] = y0.toRealRI;

    auto xhat2r = EP!(softDecision2PAM_QPSK!F, F, F)(y0r, Ar, 0.5, sigma2/2, 2);
    auto xhat2 = vector!C(N);
    xhat2[] = xhat2r.fromRealRI;

    foreach(i; 0 .. N) {
        assert(isClose(xhat1[i].re, xhat2[i].re, 1e-2));
        assert(isClose(xhat1[i].im, xhat2[i].im, 1e-2));
    }

    auto xhat3 = EP!(softDecisionQPSK!C, C, F)(y0, U, ds, V, 1, sigma2, 20);
    foreach(i; 0 .. N) {
        assert(isClose(xhat3[i].re, x0[i].re, 1e-2));
        assert(isClose(xhat3[i].im, x0[i].im, 1e-2));
    }


    auto detector = makeSVDEPDetector!(C)(QPSK!C(), U, ds, V, sigma2, 20);
    C[] xhat4;
    C[] y0arr = new C[M];
    y0arr.sliced[] = y0.sliced;

    detector.detect(y0arr, xhat4);
    foreach(i; 0 .. N) {
        assert(isClose(xhat4[i].re, x0[i].re, 1e-2));
        assert(isClose(xhat4[i].im, x0[i].im, 1e-2));
    }
}
