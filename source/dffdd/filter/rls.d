module dffdd.filter.rls;

import carbon.math;
import carbon.linear;
import std.stdio;
import std.math;
import std.complex;
import dffdd.utils.linalg;

import mir.ndslice;
//import glas.ndslice;

final class RLSAdapter(State)
{
    import std.algorithm;
    alias C = State.StateElementType;
    alias F = typeof(C.init.re);


    this(State state, real lambda, real delta = 1E-3)
    {
        _lambdaInv = 1/lambda;

        immutable size = state.state.elementsCount;
        _u = new C[size].sliced(size);
        _p = new C[size*size].sliced(size, size);
        _pu = new C[size].sliced(size);
        _uhp = new C[size].sliced(size);
        _k = new C[size].sliced(size);

        _u[] = complexZero!C;
        _p[] = complexZero!C;
        _p.diagonal[] = ((1/delta) + complexZero!C);
    }


    void adapt()(auto ref State state, C error)
    {
        immutable N = state.state.length!0;
        immutable P = state.state.length!1;

        foreach_reverse(i; P .. N*P)
            _u[i] = _u[i-P];

        _u[0 .. P] = state.state[0];

        // _pu = _p * _u
        matop_M_mul_V(_p, _u, _pu);
        _pu[] *= _lambdaInv;

        // _uhp = _u * _p
        matop_Vh_mul_M(_u, _p, _uhp);
        _uhp[] *= _lambdaInv;

        // _u * _pu
        immutable C uhpu = matop_Vh_dot_V(_u, _pu);

        _k[] = _pu[];
        _k[] /= (1 + uhpu);

        // _p = _p * _lambdaInv - _k * _uhp
        _p[] *= _lambdaInv;
        foreach(i; 0 .. _p.length!0) foreach(j; 0 .. _p.length!1)
            _p[i, j] -= _k[i] * _uhp[j];

        foreach(i; 0 .. N){
            foreach(j; 0 .. P){
                auto idx = i*P + j;
                auto kc = _k[idx].conj;
                state.weight[i, j] += kc * error;
            }
        }
    }


  private:
    real _lambdaInv;
    Slice!(C*, 1, Contiguous) _u;
    Slice!(C*, 2, Contiguous) _p;
    Slice!(C*, 1, Contiguous) _pu;
    Slice!(C*, 1, Contiguous) _uhp;
    Slice!(C*, 1, Contiguous) _k;
}


RLSAdapter!State makeRLSAdapter(State)(State state, real lambda, real delta = 1E-3)
{
    return new typeof(return)(state, lambda, delta);
}


unittest
{
   import dffdd.filter.state;

   auto state = MultiFIRState!(Complex!float)(1, 1);
   auto rlsAdapter = makeRLSAdapter(state, 0.5);
   rlsAdapter.adapt(state, Complex!float(0, 0));
}



private
void matop_M_mul_V(C)(Slice!(C*, 2, Contiguous) mat, Slice!(C*, 1, Contiguous) v, Slice!(C*, 1, Contiguous) dst)
{
    // foreach(i; 0 .. dst.length){
    //     dst[i] = complexZero!C;
    //     foreach(j; 0 .. mat.length!1)
    //         dst[i] += mat[i, j] * v[j];
    // }
    gemv(C(1, 0), mat, v, C(0, 0), dst);
}

unittest
{
    alias C = Complex!float;

    auto mat = new C[2*2].sliced(2, 2);
    auto v1 = new C[2].sliced(2);

    mat[0, 0] = C(1, 0); mat[0, 1] = C(0, 1);
    mat[1, 0] = C(2, 1); mat[1, 1] = C(2, 2);
    v1[0] = C(0, -1); v1[1] = C(1, -1);

    auto dst = new C[2].sliced(2);
    matop_M_mul_V(mat, v1, dst);
    assert(approxEqualCpx(dst[0], C(1, 0)));
    assert(approxEqualCpx(dst[1], C(5, -2)));
}


private
void matop_Vh_mul_M(C)(Slice!(C*, 1, Contiguous) v, Slice!(C*, 2, Contiguous) mat, Slice!(C*, 1, Contiguous) dst)
{
    // foreach(i; 0 .. dst.length){
    //     dst[i] = complexZero!C;
    //     foreach(j; 0 .. mat.length!0)
    //         dst[i] += v[j].conj * mat[j, i];
    // }
    gemv!"H"(C(1, 0), mat, v, C(0, 0), dst);
    foreach(i; 0 .. dst.length)
        dst[i] = dst[i].conj;
}


private
C matop_Vh_dot_V(C)(Slice!(C*, 1, Contiguous) a, Slice!(C*, 1, Contiguous) b)
{
    // C c = complexZero!C;
    // foreach(i; 0 .. a.length)
    //     c += a[i].conj * b[i];

    // return c;
    auto d = dotH(a, b);
    return d;
}
