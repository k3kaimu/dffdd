module dffdd.filter.rls;

import carbon.math;
import carbon.linear;
import std.stdio;
import std.math;
import std.complex;


final class RLSAdapter(State)
{
    import std.algorithm;

    enum bool usePower = false;

    enum size_t N = typeof(State.init.state).length;
    enum size_t P = typeof(State.init.state[0]).length;
    //static assert(P == 1);

    alias C = typeof(State.init.state[0][0]);
    alias F = typeof(C.init.re);


    this(State state, real lambda, real delta = 1E-3)
    {
        _lambdaInv = 1/lambda;

        C zero;
        zero = 0;

        _u = ones!float * zero;
        _p = identity!float * ((1/delta) + complexZero!C);
    }


    void adapt(ref State state, C error)
    {
        foreach_reverse(i; P .. N*P)
            _u[i] = _u[i-P];

        foreach(i; 0 .. P)
            _u[i] = state.state[0][i];

        auto u = _u.pref,
             uh = u.hermitian,
             p = _p.pref;

        SMatrix!(C, N*P, 1) pu = p * u * _lambdaInv;
        SMatrix!(C, 1, N*P) uhp = uh * p * _lambdaInv;
        C uhpu = uh * pu.pref;

        SMatrix!(C, N*P, 1) k = pu.pref / (1 + uhpu);
        p = p * _lambdaInv - k.pref * uhp.pref;

        //foreach(i; 0 .. N*P){
        //    auto kc = k[i].conj;
        //    state.weight[i][0] += kc * error;
        //}
        foreach(i; 0 .. N){
            foreach(j; 0 .. P){
                auto idx = i*P+j;
                auto kc = k[idx].conj;
                state.weight[i][j] += kc * error;
            }
        }
    }


  private:
    real _lambdaInv;
    SMatrix!(C, N*P, 1) _u;
    SMatrix!(C, N*P, N*P) _p;
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
