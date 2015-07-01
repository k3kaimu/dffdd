module dffdd.filter.mempoly;

import std.complex;
import std.math;

import std.stdio;

import carbon.math;
import carbon.stream;

import dffdd.filter.traits;


final class MemoryPolynomialState(C, uint N, uint Pa, uint Pb, uint Pc, uint Mb, uint Mc)
{
    import std.algorithm: max;


    alias F = typeof(C.init.re);

    enum size_t Pbc = max(Pb, Pc);
    enum size_t Pabc = max(Pa-1, Pb, Pc);
    enum size_t Mbc = Mb + Mc;
    enum size_t Mabc = Mb + Mc + 1;

    C[Mc] xmemory;
    C[Mc+1] ymemory;
    F[Pabc][Mabc] memory;
    C[Pb*Mb + Pa + Pc*Mc + N][1] weight;
    C[Pb*Mb + Pa + Pc*Mc + N][1] state;
    F[Pb*Mb + Pa + Pc*Mc + N] power;


    this()
    {
        foreach(ref e; xmemory) e = complexZero!C;
        foreach(ref e; ymemory) e = complexZero!C;
        foreach(ref e; memory) e[] = 0;
        foreach(i; 0 .. Pb*Mb + Pa + Pc*Mc + N){
            weight[0][i] = complexZero!C;
            state[0][i] = complexZero!C;
            power[i] = 0;
        }
    }


    void update(C x)
    {
        // update memory
        foreach_reverse(i; 1 .. Mabc)
            memory[i] = memory[i-1];

        auto xabs = x.abs();
        foreach(i; 0 .. Pabc)
            memory[0][i] = xabs ^^ (2*i);

        auto xn = xmemory[$-1];

        // update state[0 .. Pb*Mb + Pa + Pc*Mc]
        C[] p = state[0]; size_t idx = 0;
        foreach(i; 0 .. Mc) foreach(j; 0 .. Pc){
            p[idx] = memory[i][j] * xn;
            ++idx;
        }
        p[idx] = xn; ++idx;
        foreach(j; 1 .. Pa){
            p[idx] = memory[Mc][j-1] * xn;
            ++idx;
        }
        foreach(i; 0 .. Mb) foreach(j; 0 .. Pb){
            p[idx] = memory[Mc+1+i][j] * xn;
            ++idx;
        }

        // update state[Pb*Mb + Pa + Pc*Mc .. $]
        foreach_reverse(i; 1 .. N)
            state[0][Pb*Mb + Pa + Pc*Mc + i] = state[0][Pb*Mb + Pa + Pc*Mc + i - 1];

        foreach(i; 0 .. Pb*Mb + Pa + Pc*Mc)
            state[0][Pb*Mb + Pa + Pc*Mc] += weight[0][i] * state[0][i];

         //get max power state[Pb*Mb + Pa + Pc*Mc .. $]
        //{
        //    F maxPower = 0;
        //    foreach(i; Pb*Mb + Pa + Pc*Mc .. Pb*Mb + Pa + Pc*Mc + N)
        //    {
        //        auto c = state[0][i];
        //        auto pp = c.re^^2 + c.im^^2;
        //        if(maxPower < pp) maxPower = pp;
        //    }

        //    F maxVolt = sqrt(maxPower);
        //    //foreach(i; 0 .. Pb*Mb + Pa + Pc*Mc)
        //    //    state[0][i] /= maxVolt;

        //    foreach(i; Pb*Mb + Pa + Pc*Mc .. Pb*Mb + Pa + Pc*Mc + N)
        //    {
        //        state[0][i] /= maxVolt;
        //    }
        //}

        // update power
        foreach(i, e; state[0])
            power[i] = e.re^^2 + e.im^^2;

        // update xmemory
        foreach_reverse(i; 1 .. xmemory.length)
            xmemory[i] = xmemory[i-1];

        xmemory[0] = x;

        // update ymemory
        foreach_reverse(i; 1 .. ymemory.length)
            ymemory[i] = ymemory[i-1];
    }


    C error(C y)
    {
        ymemory[0] = y;
        C c = ymemory[Mc];
        foreach(i; Pa + Pb*Mb + Pc*Mc .. Pa + Pb*Mb + Pc*Mc + N)
            c -= weight[0][i] * state[0][i];

        //writeln(c);
        return c;
    }
}
