module dffdd.utils.linalg;

import std.complex;

nothrow: @trusted:
extern(C) int LAPACKE_zheev(
            int matrix_order,
            const char jobz, const char uplo,
            const int N, double[2]* A, const int lda,
            double* w);


extern(C) int LAPACKE_cheev(
            int matrix_order,
            const char jobz, const char uplo,
            const int N, float[2]* A, const int lda,
            float* w);



int heev(size_t N, double[2]* A, double* eigenValue)
{
    int info = LAPACKE_zheev(101, 'V', 'U', cast(uint)N, A, cast(uint)N, eigenValue);
    return info;
}


int heev(size_t N, float[2]* A, float* eigenValues)
{
    int info = LAPACKE_cheev(101, 'V', 'U', cast(uint)N, A, cast(uint)N, eigenValues);
    return info;
}
