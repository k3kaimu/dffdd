module dffdd.filter.linfir;

import std.complex;
import std.math;

import std.stdio;

import carbon.math;
import dffdd.filter.traits;


final class LinearFIR(C, size_t N)
if(N > 0)
{
    C[1][N] state;
    C[1][N] weight;


    this()
    {
        foreach(ref e; state) e = complexZero!C;
        foreach(ref e; weight) e = complexZero!C;
    }


    void update(C x)
    {
        foreach_reverse(i; 1 .. N)
            state[i] = state[i-1];

        state[0] = x;
    }


    void error(C y)
    {
        foreach(i; 0 .. N)
            y -= state[i] * weight[i];
    }
}
