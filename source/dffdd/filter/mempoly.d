module dffdd.filter.mempoly;

import std.complex;
import std.math;

import std.stdio;

import carbon.math;
import carbon.stream;

import dffdd.filter.traits;


void shiftBack(T)(T[] array) @trusted
{
    size_t n = array.length;
    auto q = array.ptr;
    auto p = q + n - 1;
    while(p != q){
        *p = *(p-1);
        --p;
    }
}


final class MemoryPolynomialState(C, size_t N, size_t P, size_t Mf, size_t Mp, bool withDCBias = true, bool withIQImbalance = true, bool usePower = true)
{
    alias F = typeof(C.re);

    C[Mf+1] ymemory;
    C[Mf+1+Mp] xmemory;
    static if(P>1) F[P-1][Mf+1+Mp] pmemory;
    C[withDCBias + (1 + (P-1)*(Mf+1+Mp)) * (withIQImbalance ? 2 : 1)][N] state;
    C[withDCBias + (1 + (P-1)*(Mf+1+Mp)) * (withIQImbalance ? 2 : 1)][N] weight;
    size_t cnt;

  static if(usePower)
  {
    F[withDCBias + (1 + (P-1)*(Mf+1+Mp)) * (withIQImbalance ? 2 : 1)] power;
  }


    this(F initP)
    {
        foreach(ref e; xmemory) e = complexZero!C;
        foreach(ref e; ymemory) e = complexZero!C;
        static if(P>1) foreach(ref e; pmemory) e = 0;
        foreach(ref e; state) foreach(ref ee; e) ee = complexZero!C;
        foreach(ref e; weight) foreach(ref ee; e) ee = complexZero!C;
        foreach(ref e; power) e = initP;
    }


    void update(C c) pure nothrow @safe @nogc
    {
        xmemory.shiftBack;
        ymemory.shiftBack;
        static if(P>1) pmemory.shiftBack;
        state.shiftBack;

        xmemory[0] = c;
        immutable cabs2 = c.re^^2 + c.im^^2;

      static if(P>1)
        foreach(p; 2 .. P+1)
            pmemory[0][p-2] = cabs2^^(p-1);

        immutable xk = xmemory[Mf];
        static if(withDCBias) state[0][0] = complexZero!C + 1;    // 1 + 0i
        state[0][withDCBias + 0] = xk;
        static if(withIQImbalance) state[0][withDCBias + 1 + (P-1)*(Mf+1+Mp)] = xk.conj;

        foreach(p; 2-2 .. P+1-2) foreach(m; 0 .. Mp+Mf+1){
          static if(P>1)
            auto xx = xk * pmemory[m][p];
          else{
            auto xx = xk;
            //assert(0);
          }

            state[0][withDCBias + 1 + p*(Mp+Mf+1) + m] = xx;

          static if(withIQImbalance)
            state[0][withDCBias + 1 + (P-1)*(Mf+1+Mp) + 1 + p*(Mf+Mp+1) + m] = xx.conj;
        }

      static if(usePower)
        foreach(i; 0 .. power.length){
            auto e = state[0][i];
            power[i] = e.re^^2 + e.im^^2;
        }
    }


    C error(C c) pure nothrow @safe @nogc
    {
        ymemory[0] = c;
        c = ymemory[Mf];
        foreach(i; 0 .. state.length) foreach(j; 0 .. state[0].length){
            c -= state[i][j] * weight[i][j];
        }

        return c;
    }


    void apply(in C[] tx, in C[] rx, C[] dst) pure nothrow @safe @nogc
    {
        foreach(ic, c; tx){
            xmemory.shiftBack;
            ymemory.shiftBack;
            static if(P>1) pmemory.shiftBack;
            state.shiftBack;

            xmemory[0] = c;
            immutable cabs2 = c.re^^2 + c.im^^2;

          static if(P>1)
            foreach(p; 2 .. P+1)
                pmemory[0][p-2] = cabs2^^(p-1);

            immutable xk = xmemory[Mf];
            static if(withDCBias) state[0][0] = complexZero!C + 1;    // 1 + 0i
            state[0][withDCBias + 0] = xk;
            static if(withIQImbalance) state[0][withDCBias + 1 + (P-1)*(Mf+1+Mp)] = xk.conj;

            foreach(p; 2-2 .. P+1-2) foreach(m; 0 .. Mp+Mf+1){
              static if(P>1)
                auto xx = xk * pmemory[m][p];
              else{
                auto xx = xk;
                //assert(0);
              }

                state[0][withDCBias + 1 + p*(Mp+Mf+1) + m] = xx;

              static if(withIQImbalance)
                state[0][withDCBias + 1 + (P-1)*(Mf+1+Mp) + 1 + p*(Mf+Mp+1) + m] = xx.conj;
            }

            ymemory[0] = rx[ic];
            dst[ic] = ymemory[Mf];
            foreach(i; 0 .. state.length) foreach(j; 0 .. state[0].length)
                dst[ic] -= state[i][j] * weight[i][j];
        }
    }
}

__EOF__


final class MemoryPolynomialStateBySIMD(C, size_t N, size_t P, size_t Mf, size_t Mp, bool withDCBias = true, bool withIQImbalance = true, bool usePower = true)
{
    alias F = typeof(C.re);

    C[Mf+1] ymemory;
    C[Mf+1+Mp] xmemory;
    //F[P-1][Mf+1+Mp] pmemory;
    Vector!(F[SIMD_N])[(Mf+1+Mp).alignSize(SIMD_N)][P-1] pmemory;
    Vector!(F[SIMD_N])[(Mf+1+Mp).alignSize(SIMD_N) * (P-1) * (withIQImbalance ? 2 : 1) + 1][N] reState, imState, reWeight, imWeight;

    size_t cnt;

  static if(usePower)
  {
    F[withDCBias + (1 + (P-1)*(Mf+1+Mp)) * (withIQImbalance ? 2 : 1)] power;
  }


    this(F initP)
    {
        foreach(ref e; xmemory) e = complexZero!C;
        foreach(ref e; ymemory) e = complexZero!C;
        foreach(ref e; pmemory) e = 0;
        foreach(ref e; reState) foreach(ref ee; e) ee[] = 0;
        foreach(ref e; imState) foreach(ref ee; e) ee[] = 0;
        foreach(ref e; reWeight) foreach(ref ee; e) ee[] = 0;
        foreach(ref e; imWeight) foreach(ref ee; e) ee[] = 0;
        static if(usePower) foreach(ref e; power) e = initP;
    }


    void update(C c) pure nothrow @safe @nogc
    {
        xmemory.shiftBack;
        ymemory.shiftBack;
        pmemory.shiftBack;
        state.shiftBack;

        xmemory[0] = c;
        immutable cabs2 = c.re^^2 + c.im^^2;

        foreach(p; 2 .. P+1)
            pmemory[0][p-2] = cabs2^^(p-1);

        immutable xk = xmemory[Mf];
        //static if(withDCBias) state[0][0] = complexZero!C + 1;    // 1 + 0i
        //state[0][withDCBias + 0] = xk;
        //static if(withIQImbalance) state[0][withDCBias + 1 + (P-1)*(Mf+1+Mp)] = xk.conj;

        //foreach(p; 2-2 .. P+1-2) foreach(m; 0 .. Mp+Mf+1){
        //    auto xx = xk * pmemory[m][p];
        //    state[0][withDCBias + 1 + p*(Mp+Mf+1) + m] = xx;

        //  static if(withIQImbalance)
        //    state[0][withDCBias + 1 + (P-1)*(Mf+1+Mp) + 1 + p*(Mf+Mp+1) + m] = xx.conj;
        //}
        foreach(p; 2-2 .. P+1-2) foreach()

      static if(usePower)
        foreach(i; 0 .. power.length){
            auto e = state[0][i];
            power[i] = e.re^^2 + e.im^^2;
        }
    }


    C error(C c) pure nothrow @safe @nogc
    {
        ymemory[0] = c;
        c = ymemory[Mf];
        foreach(i; 0 .. state.length) foreach(j; 0 .. state[0].length){
            c -= state[i][j] * weight[i][j];
        }

        return c;
    }


    void apply(in C[] tx, in C[] rx, C[] dst) pure nothrow @safe @nogc
    {
        foreach(ic, c; tx){
            xmemory.shiftBack;
            ymemory.shiftBack;
            pmemory.shiftBack;
            state.shiftBack;

            xmemory[0] = c;
            immutable cabs2 = c.re^^2 + c.im^^2;

            foreach(p; 2 .. P+1)
                pmemory[0][p-2] = cabs2^^(p-1);

            immutable xk = xmemory[Mf];
            static if(withDCBias) state[0][0] = complexZero!C + 1;    // 1 + 0i
            state[0][withDCBias + 0] = xk;
            static if(withIQImbalance) state[0][withDCBias + 1 + (P-1)*(Mf+1+Mp)] = xk.conj;

            foreach(p; 2-2 .. P+1-2) foreach(m; 0 .. Mp+Mf+1){
                auto xx = xk * pmemory[m][p];
                state[0][withDCBias + 1 + p*(Mp+Mf+1) + m] = xx;

              static if(withIQImbalance)
                state[0][withDCBias + 1 + (P-1)*(Mf+1+Mp) + 1 + p*(Mf+Mp+1) + m] = xx.conj;
            }

            ymemory[0] = rx[ic];
            dst[ic] = ymemory[Mf];
            foreach(i; 0 .. state.length) foreach(j; 0 .. state[0].length)
                dst[ic] -= state[i][j] * weight[i][j];
        }
    }
}
