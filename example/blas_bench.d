/++ dub.json:
{
    "name": "blas_bench",
    "dependencies": {
        "dffdd": { "path": ".." }
    },
    "dflags-ldc": ["-enable-cross-module-inlining", "-mcpu=native"],
}
+/
module example.blas;

import std.algorithm;
import std.datetime.stopwatch;
import std.range;
import std.stdio;

import dffdd.math.matrix;
import dffdd.math.vector;
import dffdd.math.exprtemplate;

import mir.ndslice : Slice, SliceKind;

import ldc.attributes : fastmath;

alias F = float;

extern(C) void openblas_set_num_threads(int num_threads);
extern(C) int openblas_get_num_threads();


shared static this()
{
    openblas_set_num_threads(1);
}


void main()
{
    size_t ITERMAX = 1_000_000_000;

    foreach(p; iota(1, 11)) {
        size_t M = 2^^p;
        benchmark(M, min(max(ITERMAX/M^^3, 100), 10000000));
    }

    // static foreach(M; iota(2, 16)) {
    //     benchmark(M, min(max(ITERMAX/M^^3, 100), 10000000));
    // }
}


@fastmath
void benchmark(size_t M, size_t ITER)
{
    auto A = matrix!F(M, M, 0);
    auto B = matrix!F(M, M, 0);
    auto C = matrix!F(M, M, 0);
    auto D = matrix!F(M, M, 0);
    // auto x = vector!F(M, 0);
    // auto y = vector!F(M, 0);


    A[] = 1;
    B[] = 1;
    C[] = 1;

    F a = 1.0001;
    F b = 1.0002;

    auto sA = A.sliced,
         sB = B.sliced,
         sC = C.sliced,
         sD = D.sliced;

    auto pA = A.sliced.iterator,
         pB = B.sliced.iterator,
         pC = C.sliced.iterator,
         pD = D.sliced.iterator;

    StopWatch sw;
    sw.start();
    foreach(i; 0 .. ITER) {
        //  BLAS
        // gemmBLAS(M, a, pA, pB, b, pC, pD);

        // For
        // gemmForLoopPtr(M, M, M, a, pA, pB, b, pC, pD);

        // dffdd(Lazy)
        D.noalias = a * A * B * C + b * C;

        // dffdd(Eager)
        // D.noalias = forceEvaluate(forceEvaluate(forceEvaluate(a * A) * B) + forceEvaluate(b * C));

        A[] = D;
    }
    sw.stop();
    writefln!"M = %s, %s GFLOPS"(M, M*M*(2*M+3)*1.0*ITER / sw.peek.total!"usecs" / 1e3);
}


@fastmath
void gemmForLoopPtr(F)(size_t M, size_t N, size_t K, F a, const(F)* pA, const(F)* pB, F b, const(F)* pC, F* pD)
{
    foreach(i; 0 .. M) {
        foreach(k; 0 .. K)
            foreach(j; 0 .. N)
                pD[i * N + j] += pA[i * K + k] * pB[k * N + j];
    
        foreach(j; 0 .. N)
            pD[i * N + j] = pD[i * N + j] * a + pC[i * N + j] * b;
    }
}


void gemmBLAS(size_t M, float a, const(float)* pA, const(float)* pB, float b, const(float)* pC, float* pD)
{
    import cblas;
    pD[0 .. M*M] = pC[0 .. M*M];
    cblas_sgemm(
        CBLAS_ORDER.RowMajor,
        CBLAS_TRANSPOSE.NoTrans,
        CBLAS_TRANSPOSE.NoTrans,
        cast(int)M, cast(int)M, cast(int)M,
        a,
        pA, cast(int)M,
        pB, cast(int)M,
        b,
        pD, cast(int)M);
}