module dffdd.filter.mempoly;

import std.complex;
import std.math;

import std.stdio;

import carbon.math;
import carbon.stream;

import dffdd.filter.traits;


final class MemoryPolynomialState(C, size_t N, size_t P, size_t Mf, size_t Mp, bool withDCBias = true)
{
    alias F = typeof(C.re);

    C[Mf+1] ymemory;
    C[Mf+1+Mp] xmemory;
    F[P-1][Mf+1+Mp] pmemory;
    C[withDCBias + 1 + (P-1)*(Mf+1+Mp)][N] state;
    C[withDCBias + 1 + (P-1)*(Mf+1+Mp)][N] weight;
    F[withDCBias + 1 + (P-1)*(Mf+1+Mp)] power;


    this(F initP)
    {
        foreach(ref e; xmemory) e = complexZero!C;
        foreach(ref e; ymemory) e = complexZero!C;
        foreach(ref e; pmemory) e = 0;
        foreach(ref e; state) foreach(ref ee; e) ee = complexZero!C;
        foreach(ref e; weight) foreach(ref ee; e) ee = complexZero!C;
        foreach(ref e; power) e = initP;
    }


    void update(C c)
    {
        foreach_reverse(i; 1 .. xmemory.length)
            xmemory[i] = xmemory[i-1];

        foreach_reverse(i; 1 .. ymemory.length)
            ymemory[i] = ymemory[i-1];

        foreach_reverse(i; 1 .. pmemory.length)
            pmemory[i] = pmemory[i-1];

        foreach_reverse(i; 1 .. state.length)
            state[i] = state[i-1];

        xmemory[0] = c;
        immutable cabs2 = c.re^^2 + c.im^^2;

        foreach(p; 2 .. P+1)
            pmemory[0][p-2] = cabs2^^(p-1);

        immutable xk = xmemory[Mf];
        static if(withDCBias) state[0][0] = complexZero!C + 1;    // 1 + 0i
        state[0][withDCBias ? 1 : 0] = xk;

        foreach(p; 2-2 .. P+1-2) foreach(m; 0 .. Mp+Mf+1){
            state[0][withDCBias + 1 + p*(Mp+Mf+1) + m] = xk * pmemory[m][p];
        }

        foreach(i; 0 .. power.length){
            auto e = state[0][i];
            power[i] = e.re^^2 + e.im^^2;
        }
    }


    C error(C c)
    {
        ymemory[0] = c;
        c = ymemory[Mf];
        foreach(i; 0 .. N) foreach(j; 0 .. 1 + (P-1)*(Mf+1+Mp))
            c -= state[i][j] * weight[i][j];

        return c;
    }
}

/*
final class MemoryPolynomialState(C, size_t N, size_t P, size_t Mf, size_t Mp, bool withDCBias = true)
{
    alias F = typeof(C.re);

    C[Mf+1] ymemory;
    C[Mf+1+Mp] xmemory;
    F[P-1][Mf+1+Mp] pmemory;
    C[withDCBias + 1 + (P-1)*(Mf+1+Mp)][N] state;
    C[withDCBias + 1 + (P-1)*(Mf+1+Mp)][N] weight;
    F[withDCBias + 1 + (P-1)*(Mf+1+Mp)] power;


    this(F initP)
    {
        foreach(ref e; xmemory) e = complexZero!C;
        foreach(ref e; ymemory) e = complexZero!C;
        foreach(ref e; pmemory) e = 0;
        foreach(ref e; state) foreach(ref ee; e) ee = complexZero!C;
        foreach(ref e; weight) foreach(ref ee; e) ee = complexZero!C;
        foreach(ref e; power) e = initP;
    }


    void update(C c)
    {
        foreach_reverse(i; 1 .. xmemory.length)
            xmemory[i] = xmemory[i-1];

        foreach_reverse(i; 1 .. ymemory.length)
            ymemory[i] = ymemory[i-1];

        foreach_reverse(i; 1 .. pmemory.length)
            pmemory[i] = pmemory[i-1];

        foreach_reverse(i; 1 .. state.length)
            state[i] = state[i-1];

        xmemory[0] = c;
        immutable cabs2 = c.re^^2 + c.im^^2;

        foreach(p; 2 .. P+1)
            pmemory[0][p-2] = cabs2^^(p-1);

        immutable xk = xmemory[Mf];
        static if(wid) state[0][0] = complexZero!C + 1;    // 1 + 0i
        state[0][withDCBias ? 1 : 0] = xk;

        foreach(p; 2-2 .. P+1-2) foreach(m; 0 .. Mp+Mf+1){
            state[0][withDCBias + 1 + p*(Mp+Mf+1) + m] = xk * pmemory[m][p];
        }

        foreach(i; 0 .. power.length){
            auto e = state[0][i];
            power[i] = e.re^^2 + e.im^^2;
        }
    }


    C error(C c)
    {
        ymemory[0] = c;
        c = ymemory[Mf];
        foreach(i; 0 .. N) foreach(j; 0 .. 1 + (P-1)*(Mf+1+Mp))
            c -= state[i][j] * weight[i][j];

        return c;
    }
}
*/