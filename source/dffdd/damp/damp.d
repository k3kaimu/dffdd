module dffdd.damp.damp;

import std.math : sqrt, isClose, SQRT2, SQRT1_2;
import std.numeric : stddot = dotProduct;
import std.algorithm : stdmap = map;
import std.traits : isFloatingPoint, isArray;
import std.typecons : Tuple, tuple;

import dffdd.math.vector;
import dffdd.math.matrix;
import dffdd.math.exprtemplate : hasMemoryView;

import mir.ndslice : SliceKind, Contiguous;


Tuple!(F, "value", F, "drv") softThrBayesOpt(alias arrP_, alias arrR_, F)(F x, F sigma, F delta)
if(isArray!(typeof(arrP_)) && isArray!(typeof(arrR_)) && isFloatingPoint!F && arrP_.length >= 2 && arrP_.length == arrR_.length)
{
    import std.math : isNaN, abs;
    import std.algorithm : min;
    import dffdd.math.math : fast_exp;

    enum F[] arrP = cast(F[])arrP_;
    enum F[] arrR = arrR_;

    version(none)
    {
        F xr2_min = F.infinity;
        static foreach(i; 0 .. arrP.length) {
            xr2_min = min((arrR[i] - x)^^2, xr2_min);
        }
    }
    else
    {
        enum F xr2_min = F(0);
    }

    immutable expCoef = -0.5 * delta / sigma^^2;

    F p_e = F(0),
      p_r_e = F(0),
      p_r2_e = F(0);
    static foreach(i; 0 .. arrP.length) {{
        F expvalue = fast_exp!F(expCoef * ((arrR[i] - x)^^2 - xr2_min));
        if(expvalue.isNaN) expvalue = 1;

        p_e += arrP[i] * expvalue;
        p_r_e += arrP[i] * arrR[i] * expvalue;
        p_r2_e += arrP[i] * arrR[i] * arrR[i] * expvalue;
    }}

    immutable drv = -2 * expCoef * (p_r2_e * p_e - p_r_e^^2) / p_e^^2;
    return typeof(return)(
        p_r_e / p_e,
        (abs(sigma) > 1e-15 && !drv.isNaN) ? drv : 0
    );
}

unittest
{
    assert(softThrBayesOpt!([0.5, 0, 0.5], [-1, 0, 1])(0.001, 0.1, 0.8).value.isClose(0.0798298, 1e-4));
    assert(softThrBayesOpt!([0.5, 0, 0.5], [-1, 0, 1])(0.001, 0.1, 0.8).drv.isClose(79.4902, 1e-4));
}


alias softThrBayesOpt2PAM_BPSK(F) = softThrBayesOpt!([0.5, 0.5], [-1, 1], F);
alias softThrBayesOpt2PAM_QPSK(F) = softThrBayesOpt!([0.5, 0.5], [-SQRT1_2, SQRT1_2], F);
alias softThrBayesOpt4PAM_16QAM(F) = softThrBayesOpt!([0.25, 0.25, 0.25, 0.25], [-3/sqrt(10.0L), -1/sqrt(10.0L), 1/sqrt(10.0L), 3/sqrt(10.0L)], F);
alias softThrBayesOpt8PAM_64QAM(F) = softThrBayesOpt!([0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125], [-7/sqrt(42.0L), -5/sqrt(42.0L), -3/sqrt(42.0L), -1/sqrt(42.0L), 1/sqrt(42.0L), 3/sqrt(42.0L), 5/sqrt(42.0L), 7/sqrt(42.0L)], F);

pragma(msg, typeof(softThrBayesOpt8PAM_64QAM!float));

// Tuple!(F[], "x_hat", F[], "arrMSE")

Vector!(F, Contiguous)
    BODAMP(alias prox, F, VecY, MatA)
        (in VecY y_, in MatA A_, in F delta, in F theta0, size_t nIteration)
if(isFloatingPoint!F && isVectorLike!VecY && hasMemoryView!VecY && isMatrixLike!MatA && hasMemoryView!MatA && is(VecY.ElementType : F) && is(MatA.ElementType : F))
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

    alias ProxResult = typeof(prox(F(0), theta0, delta));
    auto wprox = slice!ProxResult(N);
    wprox[] = w.sliced.map!(e => prox(e, theta, delta));

    foreach(iter; 1 .. nIteration) {
        z[] = y - A * x + (1/delta) * z * wprox.map!"a.drv".mean;
        theta = sqrt(norm2(z) / N);
        w[] = x + A.T * z;
        wprox[] = w.sliced.map!(e => prox(e, theta, delta));
        x[] = wprox.map!"a.value".vectored;
    }

    return x;
}


unittest
{
    import std.stdio;


    import mir.ndslice;

    auto A = [0.154508497187474,  -0.404508497187474,   0.500000000000000, -0.404508497187474,  0.154508497187474,    0.154508497187474,  -0.404508497187474,   0.500000000000000, -0.404508497187473,  0.154508497187473,
    0.226995249869773,  -0.493844170297569,   0.353553390593274, 0.0782172325201156,  -0.445503262094184,   0.445503262094184, -0.0782172325201147, -0.353553390593274,  0.493844170297569,    -0.226995249869773,
    0.445503262094184,  0.0782172325201155,   -0.353553390593274,    -0.493844170297569, -0.226995249869773,  0.226995249869773,    0.493844170297569,  0.353553390593274,    -0.0782172325201148,    -0.445503262094184,
    0.475528258147577,  0.293892626146237,    3.06161699786838e-17,   -0.293892626146237,    -0.475528258147577, -0.475528258147577,  -0.293892626146237,   -9.18485099360515e-17, 0.293892626146236,   0.475528258147577,
    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,
    0.293892626146237,  -0.475528258147577,   -9.18485099360515e-17, 0.475528258147577,   -0.293892626146236,    -0.293892626146237, 0.475528258147577,   2.75545529808155e-16,  -0.475528258147577,   0.293892626146237,
    0.493844170297569,  0.445503262094184,    0.353553390593274,  0.226995249869773,    0.0782172325201155, -0.0782172325201154, -0.226995249869773,  -0.353553390593274,   -0.445503262094184,    -0.493844170297569,
    0.353553390593274,  -0.353553390593274,   -0.353553390593274,    0.353553390593274,  0.353553390593274,    -0.353553390593273, -0.353553390593274,  0.353553390593273,    0.353553390593274,  -0.353553390593273,
    ].sliced(8, 10).matrixed;

    auto y = [
        -1.81199936706226,
        -2.35139170510204,
        -1.10222286082885,
        -0.555771890670480,
        0.0283159355873546,
        0.911307827117349,
        0.0769074757814163,
        -0.686546300398577,].sliced.vectored;

    auto delta = 0.8;
    enum arrP = [0.5, 0, 0.5];
    enum arrR = [-1, 0, 1];
    auto theta0 = sqrt(stddot(arrP, arrR.stdmap!"a^^2"));

    int nIteration = 2;
    import std.stdio;

    auto xhat1 = BODAMP!((e, v, d) => softThrBayesOpt!(arrP, arrR)(e, v, d))(y, A, delta, theta0, 3);
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

    auto xhat2 = BODAMP!((e, v, d) => softThrBayesOpt!(arrP, arrR)(e, v, d))(y, A, delta, theta0, 20);
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
}


// unittest
// {
//     import std.stdio;


//     import mir.ndslice;

//     auto A = [0.154508497187474,  -0.404508497187474,   0.500000000000000, -0.404508497187474,  0.154508497187474,    0.154508497187474,  -0.404508497187474,   0.500000000000000, -0.404508497187473,  0.154508497187473,
//     0.226995249869773,  -0.493844170297569,   0.353553390593274, 0.0782172325201156,  -0.445503262094184,   0.445503262094184, -0.0782172325201147, -0.353553390593274,  0.493844170297569,    -0.226995249869773,
//     0.445503262094184,  0.0782172325201155,   -0.353553390593274,    -0.493844170297569, -0.226995249869773,  0.226995249869773,    0.493844170297569,  0.353553390593274,    -0.0782172325201148,    -0.445503262094184,
//     0.475528258147577,  0.293892626146237,    3.06161699786838e-17,   -0.293892626146237,    -0.475528258147577, -0.475528258147577,  -0.293892626146237,   -9.18485099360515e-17, 0.293892626146236,   0.475528258147577,
//     0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,    0.353553390593274,  0.353553390593274,
//     0.293892626146237,  -0.475528258147577,   -9.18485099360515e-17, 0.475528258147577,   -0.293892626146236,    -0.293892626146237, 0.475528258147577,   2.75545529808155e-16,  -0.475528258147577,   0.293892626146237,
//     0.493844170297569,  0.445503262094184,    0.353553390593274,  0.226995249869773,    0.0782172325201155, -0.0782172325201154, -0.226995249869773,  -0.353553390593274,   -0.445503262094184,    -0.493844170297569,
//     0.353553390593274,  -0.353553390593274,   -0.353553390593274,    0.353553390593274,  0.353553390593274,    -0.353553390593273, -0.353553390593274,  0.353553390593273,    0.353553390593274,  -0.353553390593273,
//     ].sliced(8, 10).matrixed;

//     auto y = [
//         -1.28215055669155,
//         -1.66816398020785,
//         -0.754100947195077,
//         -0.383613576155461,
//         0.0283159355873546,
//         0.632749822786152,
//         0.0594063301905325,
//         -0.479439519212030].sliced.vectored;

//     auto delta = 0.8;
//     enum arrP = [0.5, 0, 0.5];
//     enum double[] arrR = [-SQRT1_2, 0, SQRT1_2];
//     auto theta0 = sqrt(stddot(arrP, arrR.stdmap!"a^^2"));

//     int nIteration = 2;
//     import std.stdio;

//     auto xhat1 = BODAMP!((e, v, d) => softThrBayesOpt!(arrP, arrR)(e, v, d))(y, A, delta, theta0, 3);
//     writeln(xhat1.sliced);
//     // assert(xhat1[0].isClose(-0.9729, 1e-3));
//     // assert(xhat1[1].isClose(+0.9744, 1e-3));
//     // assert(xhat1[2].isClose(-0.9220, 1e-3));
//     // assert(xhat1[3].isClose(+0.9716, 1e-3));
//     // assert(xhat1[4].isClose(+0.8439, 1e-3));
//     // assert(xhat1[5].isClose(-0.9594, 1e-3));
//     // assert(xhat1[6].isClose(+0.9440, 1e-3));
//     // assert(xhat1[7].isClose(-0.8160, 1e-3));
//     // assert(xhat1[8].isClose(-0.9436, 1e-3));
//     // assert(xhat1[9].isClose(+0.8961, 1e-3));

//     // auto xhat2 = BODAMP!((e, v, d) => softThrBayesOpt!(arrP, arrR)(e, v, d))(y, A, delta, theta0, 20);
//     // assert(xhat2[0].isClose(-1, 1e-3));
//     // assert(xhat2[1].isClose(+1, 1e-3));
//     // assert(xhat2[2].isClose(-1, 1e-3));
//     // assert(xhat2[3].isClose(+1, 1e-3));
//     // assert(xhat2[4].isClose(+1, 1e-3));
//     // assert(xhat2[5].isClose(-1, 1e-3));
//     // assert(xhat2[6].isClose(+1, 1e-3));
//     // assert(xhat2[7].isClose(-1, 1e-3));
//     // assert(xhat2[8].isClose(-1, 1e-3));
//     // assert(xhat2[9].isClose(+1, 1e-3));
// }
