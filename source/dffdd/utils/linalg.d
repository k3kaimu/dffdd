module dffdd.utils.linalg;

import std.experimental.ndslice;
import std.complex;


enum Order
{
    RowMajor = 101,
    ColMajor = 102,
}

nothrow @trusted extern(C)
{
    int LAPACKE_zheev(
            int matrix_order,
            const char jobz, const char uplo,
            const int N, double[2]* A, const int lda,
            double* w);


    int LAPACKE_cheev(
            int matrix_order,
            const char jobz, const char uplo,
            const int N, float[2]* A, const int lda,
            float* w);


    int LAPACKE_cgelss(
            int matrix_order, int m, int n,
            int nrhs, float[2]* a,
            int lda, float[2]* b,
            int ldb, float* s, float rcond,
            int* rank );

    int LAPACKE_zgelss(
            int matrix_order, int m, int n,
            int nrhs, double[2]* a,
            int lda, double[2]* b,
            int ldb, double* s, double rcond,
            int* rank );

    int LAPACKE_cgetrf(
            int matrix_order, int m, int n,
            float[2]* a, int lda, int* ipiv
        );

    int LAPACKE_zgetrf(
            int matrix_order, int m, int n,
            double[2]* a, int lda, int* ipiv
        );

    int LAPACKE_cgetri(
            int matrix_order, int n,
            float[2]* a, int lda, const int* ipiv
        );

    int LAPACKE_zgetri(
            int matrix_order, int n,
            double[2]* a, int lda, const int* ipiv
        );
}


private
bool checkMatrixStride(E)(in ref Slice!(2, E*) mat)
{
    return mat.stride!0 == 1 || mat.stride!1 == 1;
}


unittest
{
    auto mat33 = new int[9].sliced(3, 3);
    assert(checkMatrixStride(mat33));

    auto mat22 = mat33[0 .. 2, 1 .. 3];
    assert(checkMatrixStride(mat22));
}


private
bool isTransposed(E)(in ref Slice!(2, E*) mat)
{
    return mat.stride!0 == 1;
}

unittest
{
    auto mat33 = new int[9].sliced(3, 3);
    assert(!mat33.isTransposed);

    mat33 = mat33.transposed;
}


private
size_t leadingDimension(E)(in ref Slice!(2, E*) mat)
{
    if(mat.isTransposed)
        return mat.stride!1;
    else
        return mat.stride!0;
}


/*
int heev(size_t N, double[2]* A, double* eigenValue)
{
    int info = LAPACKE_zheev(Order.RowMajor, 'V', 'U', cast(uint)N, A, cast(uint)N, eigenValue);
    return info;
}


int heev(size_t N, float[2]* A, float* eigenValues)
{
    int info = LAPACKE_cheev(Order.RowMajor, 'V', 'U', cast(uint)N, A, cast(uint)N, eigenValues);
    return info;
}*/


int heev(C, R)(Slice!(2, C*) matA, R[] ev)
if(is(typeof(C.init.re) == R))
in{
    assert(matA.checkMatrixStride);
    assert(!matA.isTransposed);
    assert(ev.length == matA.length!0);
}
body{
  static if(is(R == float))
    return LAPACKE_cheev(Order.RowMajor, 'V', 'U', cast(uint)matA.length!0, cast(R[2]*)matA.ptr, matA.leadingDimension, ev.ptr);
  else
    return LAPACKE_zheev(Order.RowMajor, 'V', 'U', cast(uint)matA.length!0, cast(R[2]*)matA.ptr, matA.leadingDimension, ev.ptr);
}


int gelss(C, R)(Slice!(2, C*) matA, Slice!(1, C*) b, R[] sv, int* rank = null)
if(is(typeof(C.init.re) == R))
in{
    assert(matA.checkMatrixStride);
    assert(!matA.isTransposed);
    assert(sv.length >= min(matA.length!0, matA.length!1));
}
body{
    int rankDest;

  static if(is(R == float))
    int info = LAPACKE_cgelss(Order.RowMajor, matA.length!0, matA.length!1, 1, cast(R[2]*)matA.ptr, matA.leadingDimension, cast(R[2]*)b.ptr, 1, sv.ptr, -1, &rankDest);
  else
    int info = LAPACKE_zgelss(Order.RowMajor, matA.length!0, matA.length!1, 1, cast(R[2]*)matA.ptr, matA.leadingDimension, cast(R[2]*)b.ptr, 1, sv.ptr, -1, &rankDest);

    if(rank)
        *rank = rankDest;

    return info;
}


E dot(E)(Slice!(1, E*) x, Slice!(1, E*) y)
in{
    assert(x.length == y.length);
}
body{
  static if(is(E == float))
    return cblas_sdot(x.length, x.ptr, x.stride!0, y.ptr, y.stride!0);
  else static if(is(E == double))
    return cblas_ddot(x.length, x.ptr, x.stride!0, y.ptr, y.stride!0);
  else static if(is(typeof(E.init.re) == float))
  {
    auto ret = cblas_cdotu(x.length, x.ptr, cast(_cfloat)x.stride!0, y.ptr, cast(_cfloat)y.stride!0);

    static if(is(E == cfloat))
        return ret.re + ret.im*1.0fi;
    else
        return E(ret.re, ret.im);
  }
  else static if(is(typeof(E.init.re) == double))
  {
    auto ret = cblas_cdotu(x.length, x.ptr, cast(_cdouble)x.stride!0, y.ptr, cast(_cdouble)y.stride!0);

    static if(is(E == cfloat))
        return ret.re + ret.im*1.0i;
    else
        return E(ret.re, ret.im);
  }
  else
    static assert(0);
}


alias dotT = dot;

E dotH(E)(Slice!(1, E*) x, Slice!(1, E*) y)
in{
    assert(x.length == y.length);
}
body{
  static if(is(E == float) || is(E == double))
    return dot(x, y);
  else static if(is(typeof(E.init.re) == float))
  {
    auto ret = cblas_cdotc(x.length, x.ptr, cast(_cfloat)x.stride!0, y.ptr, cast(_cfloat)y.stride!0);

    static if(is(E == cfloat))
        return ret.re + ret.im*1.0fi;
    else
        return E(ret.re, ret.im);
  }
  else static if(is(typeof(E.init.re) == double))
  {
    auto ret = cblas_cdotc(x.length, x.ptr, cast(_cdouble)x.stride!0, y.ptr, cast(_cdouble)y.stride!0);

    static if(is(E == cfloat))
        return ret.re + ret.im*1.0i;
    else
        return E(ret.re, ret.im);
  }
  else
    static assert(0);
}


template transposeValueImpl(string op)
if(op == "" || op == "*" || op == "T" || op == "H")
{
  static if(op == "")
    enum transposeValue = Transpose.NoTrans;
  else static if(op == "*")
    enum transposeValue = Transpose.ConjNoTrans;
  else static if(op == "T")
    enum transposeValue = Transpose.Trans;
  else static if(op == "H")
    enum transposeValue = Transpose.ConjTrans;
}


Transpose transposeValue(string op)(Transpose value, bool isTransposed)
{
    if(!isTransposed)
        return transposeValueImpl!op;

    final switch(transposeValueImpl!op)
    {
      case Transpose.NoTrans:
        return Transpose.Trans;
      case Transpose.Trans:
        return Transpose.NoTrans;
      case Transpose.ConjNoTrans:
        return Transpose.ConjTrans;
      case Transpose.ConjTrans:
        return Transpose.ConjNoTrans;
    }
}


void gemv(string opA = "", E)(E alpha, Slice!(2, E*) matA, Slice!(1, E*) x, E beta, Slice!(1, E*) y)
{
    immutable isTA = matA.isTransposed;

    if(isTA)
        matA = matA.trasnposed;

  static if(is(E == float))
    cblas_sgemv(Order.RowMajor, transposeValue!opA(isTA), matA.length!0, matA.length!1, alpha, matA.ptr, matA.stride!0, x.ptr, x.stride!0, beta, y.ptr, y.stride!0);
  else static if(is(E == double))
    cblas_dgemv(Order.RowMajor, transposeValue!opA(isTA), matA.length!0, matA.length!1, alpha, matA.ptr, matA.stride!0, x.ptr, x.stride!0, beta, y.ptr, y.stride!0);
  else static if(is(typeof(E.init.re) == float))
    cblas_cgemv(Order.RowMajor, transposeValue!opA(isTA), matA.length!0, matA.length!1, alpha, cast(_cfloat*)matA.ptr, matA.stride!0, cast(_cfloat*)x.ptr, x.stride!0, beta, cast(_cfloat*)y.ptr, y.stride!0);
  else static if(is(typeof(E.init.re) == double))
    cblas_zgemv(Order.RowMajor, transposeValue!opA(isTA), matA.length!0, matA.length!1, alpha, cast(_cdouble*)matA.ptr, matA.stride!0, cast(_cdouble*)x.ptr, x.stride!0, beta, cast(_cdouble*)y.ptr, y.stride!0);
  else
    static assert(0);
}


/**
逆行列
*/
void gemi(E)(Slice!(2, E*) matA)
in{
    assert(!matA.isTransposed);
    assert(matA.length!0 == matA.length!1);
}
body{
    auto ipiv = new int[matA.length!0];

    immutable n = matA.length!0;
    immutable lda = matA.stride!0;

  static if(is(typeof(E.init.re) == float))
  {
    LAPACKE_cgetrf(Order.RowMajor, cast(uint)n, cast(uint)n, cast(float[2]*)matA.ptr, cast(uint)lda, ipiv.ptr);
    LAPACKE_cgetri(Order.RowMajor, cast(uint)n, cast(float[2]*)matA.ptr, cast(uint)lda, ipiv.ptr);
  }
  else static if(is(typeof(E.init.re) == double))
  {
    LAPACKE_zgetrf(Order.RowMajor, cast(uint)n, cast(uint)n, cast(double[2]*)matA.ptr, cast(uint)lda, ipiv.ptr);
    LAPACKE_zgetri(Order.RowMajor, cast(uint)n, cast(double[2]*)matA.ptr, cast(uint)lda, ipiv.ptr);
  }
  else static assert(0);
}

unittest
{
    import std.math;
    import std.stdio;
    auto mat = new Complex!float[4].sliced(2, 2);
    mat[0, 0] = 1;
    mat[0, 1] = 0;
    mat[1, 0] = 0;
    mat[1, 1] = 1;

    gemi(mat);
    assert(approxEqual(mat[0, 0].re, 1));
    assert(approxEqual(mat[0, 1].re, 0));
    assert(approxEqual(mat[1, 0].re, 0));
    assert(approxEqual(mat[1, 1].re, 1));
    foreach(i; [0, 1]) foreach(j; [0, 1]) assert(approxEqual(mat[i, j].im, 0));


    mat[0, 0] = 3;
    mat[0, 1] = 4;
    mat[1, 0] = 1;
    mat[1, 1] = 2;
    gemi(mat);
    //writeln(mat);

    assert(approxEqual(mat[0, 0].re, 1));
    assert(approxEqual(mat[0, 1].re, -2));
    assert(approxEqual(mat[1, 0].re, -0.5));
    assert(approxEqual(mat[1, 1].re, 1.5));
    foreach(i; [0, 1]) foreach(j; [0, 1]) assert(approxEqual(mat[i, j].im, 0));
}


/**
y[j] = sum_i mx[i, j] * a[i] のa[i]を最小二乗法で求める．
結果は，y[0 .. P]に上書きされる．(P: mx.length!0)
*/
Complex!float[] leastSquareEstimate(Slice!(2, Complex!float*) mx, Complex!float[] y)
{
    import std.algorithm : min, max;

    static float[] sworkSpace;

    immutable L = y.length;
    immutable P = mx.length!0;

    int rankN = void;
    static void adjustSize(T)(ref T[] arr, size_t n) { if(arr.length < n) arr.length = n; }

    adjustSize(sworkSpace, min(L, P));
    LAPACKE_cgelss(102, cast(int)L, cast(int)P, 1, cast(float[2]*)&(mx[0, 0]), cast(int)L, cast(float[2]*)y.ptr, cast(int)max(L, P), sworkSpace.ptr, 0.00001f, &rankN);

    return y[0 .. P];
}


__EOF__

/* 
 * svdcomp - SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs_and_types.h"
 
static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


int dsvd(float **a, int m, int n, float *w, float **v)
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs((double)a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] = (float)((double)a[k][i]/scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                f = (double)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (float)(f - g);
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += (float)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] = (float)((double)a[k][i]*scale);
            }
        }
        w[i] = (float)(scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs((double)a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] = (float)((double)a[i][k]/scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                f = (double)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (float)(f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++) 
                            a[j][k] += (float)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = (float)((double)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j][i] = (float)(((double)a[i][j] / (double)a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++) 
                        v[k][j] += (float)(s * (double)v[k][i]);
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += (float)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = (float)((double)a[j][i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (float)h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (float)(y * c + z * s);
                            a[j][i] = (float)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = (float)(-z);
                    for (j = 0; j < n; j++) 
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (float)(x * c + z * s);
                    v[jj][i] = (float)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (float)z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (float)(y * c + z * s);
                    a[jj][i] = (float)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (float)x;
        }
    }
    free((void*) rv1);
    return(1);
}
