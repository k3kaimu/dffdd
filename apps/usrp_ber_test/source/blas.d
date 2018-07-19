module blas;

import std.complex;

//import cblas;

extern(C)
{
    alias CBLAS_INDEX = size_t;
    alias blasint = int;

    CBLAS_INDEX cblas_isamax(in blasint n, in float  *x, in blasint incx);
    CBLAS_INDEX cblas_idamax(in blasint n, in double *x, in blasint incx);
    CBLAS_INDEX cblas_icamax(in blasint n, in float  *x, in blasint incx);
    CBLAS_INDEX cblas_izamax(in blasint n, in double *x, in blasint incx);
}


template BLASImpl()
{
    CBLAS_INDEX ixamax(size_t n, in float* x, size_t incx = 1) { return cblas_isamax(cast(blasint)n, x, cast(blasint) incx); }
    CBLAS_INDEX ixamax(size_t n, in double* x, size_t incx = 1) { return cblas_idamax(cast(blasint)n, x, cast(blasint) incx); }
    CBLAS_INDEX ixamax(size_t n, in Complex!float* x, size_t incx = 1) { return cblas_icamax(cast(blasint)n, cast(float*)x, cast(blasint) incx); }
    CBLAS_INDEX ixamax(size_t n, in Complex!double* x, size_t incx = 1) { return cblas_izamax(cast(blasint)n, cast(double*)x, cast(blasint) incx); }

    CBLAS_INDEX ixamax(in float[] x) { return cblas_isamax(cast(blasint)x.length, x.ptr, 1); }
    CBLAS_INDEX ixamax(in double[] x) { return cblas_idamax(cast(blasint)x.length, x.ptr, 1); }
    CBLAS_INDEX ixamax(in Complex!float[] x) { return cblas_icamax(cast(blasint)x.length, cast(float*)x.ptr, 1); }
    CBLAS_INDEX ixamax(in Complex!double[] x) { return cblas_izamax(cast(blasint)x.length, cast(double*)x.ptr, 1); }
}

alias BLAS = BLASImpl!();
