module dffdd.filter.diagonal;

import std.complex;
import std.math;

import std.stdio;

import carbon.stream;
import carbon.traits;

import dffdd.filter.traits;


final class DiagonalState(C, alias polyTerm, alias polyTermABS, uint P, uint N)
{
    alias F = typeof(C.init.re);

    C[P][N] weight;
    C[P][N] state;
    F[P] power;


    this()
    {
        C zero;
      static if(isBuiltInComplex!C)
        zero = 0+0i;
      else
        zero = C(0, 0);

        foreach(ref e; weight) foreach(ref ee; e) ee = zero;
        foreach(ref e; state) foreach(ref ee; e) ee = zero;
        foreach(ref e; power) e = 0;
    }


    void update(C c)
    {
        foreach_reverse(i, ref e; state)
            if(i != 0) e = state[i-1];

        auto cabs = c.abs();
        foreach(p, ref e; state[0])
            e = polyTerm(c, cabs, p);

        foreach(p, ref e; power)
            e = polyTermABS(cabs, p);
    }


    C error(C y)
    {
        C c = y;

        foreach(n; 0 .. N) foreach(p; 0 .. P)
            c -= weight[n][p] * state[n][p];

        return c;
    }
}
