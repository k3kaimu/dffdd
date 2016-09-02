module dffdd.filter.rls;

import carbon.math;
import carbon.linear;
import std.stdio;
import std.math;
import std.complex;
import std.experimental.ndslice;


final class RLSAdapter(State)
{
    import std.algorithm;

    enum bool usePower = false;

    //enum size_t N = typeof(State.init.state).length;
    //enum size_t P = typeof(State.init.state[0]).length;
    //static assert(P == 1);

    //alias C = typeof(State.init.state[0][0]);
    alias C = State.StateElementType;
    alias F = typeof(C.init.re);


    this(State state, real lambda, real delta = 1E-3)
    {
        _lambdaInv = 1/lambda;

        C zero = complexZero!C;

        immutable size = state.state.elementsCount;

        // allocate
        /*
        _u = matrix!(C, dynamic, 1)(size);
        _p = matrix!(C, dynamic, dynamic)(size, size);
        _pu = matrix!(C, dynamic, 1)(size);
        _uhp = matrix!(C, 1, dynamic)(size);
        _k = matrix!(C, dynamic, 1)(size);

        // initialize
        _u = ones!float * zero;
        _p = identity!float * ((1/delta) + complexZero!C);
        */
        _u = new C[size].sliced(size);
        _p = new C[size*size].sliced(size, size);
        _pu = new C[size].sliced(size);
        _uhp = new C[size].sliced(size);
        _k = new C[size].sliced(size);

        _u[] = complexZero!C;
        _p[] = complexZero!C;
        _p.diagonal[] = ((1/delta) + complexZero!C);
    }


    void adapt(ref State state, C error)
    {
        immutable N = state.state.length!0;

      static if(state.state.shape.length > 1)
        immutable P = state.state.length!1;
      else 
        immutable P = 1;

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
        C uhpu = matop_Vh_dot_V(_u, _pu);

        _k[] = _pu[];
        _k[] /= (1 + uhpu);

        // _p = _p * _lambdaInv - _k * _uhp
        _p[] *= _lambdaInv;
        foreach(i; 0 .. _p.length!0) foreach(j; 0 .. _p.length!1)
            _p[i, j] -= _k[i] * _uhp[j];

        foreach(i; 0 .. N){
          static if(state.state.shape.length > 1)
          {
            foreach(j; 0 .. P){
                auto idx = i*P + j;
                auto kc = _k[idx].conj;
                state.weight[i, j] += kc * error;
            }
          }
          else
          {
                auto idx = i;
                auto kc = _k[idx].conj;
                state.weight[i] += kc * error;
          }
        }
    }


  private:
    real _lambdaInv;
    Slice!(1, C*) _u;
    Slice!(2, C*) _p;
    Slice!(1, C*) _pu;
    Slice!(1, C*) _uhp;
    Slice!(1, C*) _k;
}


RLSAdapter!State makeRLSAdapter(State)(State state, real lambda, real delta = 1E-3)
{
    return new typeof(return)(state, lambda, delta);
}


unittest
{
    import dffdd.filter.mempoly;

    auto state = new PowerState!(cfloat, 8, 1)(1);
    auto rlsAdapter = makeRLSAdapter(state, 0.5);
}


private
void matop_M_mul_V(C)(Slice!(2, C*) mat, Slice!(1, C*) v, Slice!(1, C*) dst)
{
    foreach(i; 0 .. dst.length){
        dst[i] = complexZero!C;
        foreach(j; 0 .. mat.length!1)
            dst[i] += mat[i, j] * v[j];
    }
}


private
void matop_Vh_mul_M(C)(Slice!(1, C*) v, Slice!(2, C*) mat, Slice!(1, C*) dst)
{
    foreach(i; 0 .. dst.length){
        dst[i] = complexZero!C;
        foreach(j; 0 .. mat.length!0)
            dst[i] += v[j].conj * mat[j, i];
    }
}


private
C matop_Vh_dot_V(C)(Slice!(1, C*) a, Slice!(1, C*) b)
{
    C c = complexZero!C;
    foreach(i; 0 .. a.length)
        c += a[i].conj * b[i];

    return c;
}