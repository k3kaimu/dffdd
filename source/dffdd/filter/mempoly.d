module dffdd.filter.mempoly;

import std.complex;
import std.math;
import std.typecons;

import std.stdio;

import carbon.math;
import carbon.stream;

import dffdd.filter.traits;


enum FilterOptions
{
    none,
    withDCBias       = 1 << 0,
    withConjugate    = 1 << 1,
    useConjugate     = 1 << 2,
    usePower         = 1 << 3,
}


bool hasOption(string name)(ubyte value)
{
    return mixin("(value & FilterOptions." ~ name ~ ") != 0");
}

//alias FilterSpec = std.typecons.BitFlags!FilterOptions;

//alias WithDCBias = Flag!"withDCBias";
//alias WithConjugate = Flag!"withConjugate";
//alias UseConjugateOnly = Flag!"useConjugateOnly";
//alias UsePower = Flag!"usePower";


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


final class FIRFilter(C, size_t N, bool usePower = false)
if(N > 0)
{
    C[1][N] state;
    C[1][N] weight;
    static if(usePower) typeof(C.re)[1] power;


    this()
    {
        foreach(ref e; state) e = complexZero!C;
        foreach(ref e; weight) e = complexZero!C;
        static if(usePower) power[0] = 1;
    }


    void update(C x)
    {
        foreach_reverse(i; 1 .. N)
            state[i] = state[i-1];

        state[0] = x;
        static if(usePower) power[0] = x.re ^^2 + x.im ^^ 2;
    }


    C error(C y)
    {
        foreach(i; 0 .. N)
            y -= state[i][0] * weight[i][0];

        return y;
    }
}


final class MemoryPolynomialState(C, size_t N, size_t P, size_t Mf, size_t Mp, bool withDCBias = true, bool withIQImbalance = true, bool usePower = true, size_t startP = 0)
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
        foreach(i; 0 .. state.length) foreach(j; startP .. state[0].length){
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


final class MultiBasisFunctionsState(C, size_t N, bool usePower, Funcs...)
if(Funcs.length > 0)
{
    alias F = typeof(C.re);
    enum size_t numOfFuncs = Funcs.length;

    C[numOfFuncs][N] state;
    C[numOfFuncs][N] weight;
    size_t cnt;

  static if(usePower)
  {
    F[numOfFuncs] power;
  }


    this(F initP)
    {
        //foreach(ref e; xmemory) e = complexZero!C;
        //foreach(ref e; ymemory) e = complexZero!C;
        foreach(ref e; state) foreach(ref ee; e) ee = complexZero!C;
        foreach(ref e; weight) foreach(ref ee; e) ee = complexZero!C;
        static if(usePower) foreach(ref e; power) e = initP;
    }


    void update(C xk) pure nothrow @safe @nogc
    {
        state.shiftBack;

        foreach(fi, f; Funcs)
            state[0][fi] = f(xk);

      static if(usePower)
        foreach(i, ref pe; power){
            auto e = state[0][i];
            pe = e.re^^2 + e.im^^2;
        }
    }


    C error(C c) pure nothrow @safe @nogc
    {
        foreach(i; 0 .. N) foreach(j; 0 .. numOfFuncs)
            c -= state[i][j] * weight[i][j];

        return c;
    }


    void apply(in C[] tx, in C[] rx, C[] dst) pure nothrow @safe @nogc
    {
        foreach(ic, c; tx){
            update(c);
            dst[ic] = error(rx[ic]);
        }
    }
}


template PowerState(C, size_t N, size_t P, ubyte filterSpec = 0)
if(!filterSpec.hasOption!"withDCBias")
{
    import std.format;
    import std.array;

    string buildArgs()
    {
        auto app = appender!string();

        foreach(i; 0 .. (P+1)/2){
            immutable p = filterSpec.hasOption!"useConjugate" ? "x.conj" : "x";

            app.formattedWrite("(x => %s * (x.re^^2 + x.im^^2)^^%s), ", p, i*2);

            if(filterSpec.hasOption!"withConjugate")
                app.formattedWrite("(x => %s.conj * (x.re^^2 + x.im^^2)^^%s), ", p, i*2);
        }

        return app.data;
    }


    mixin(format("alias PowerState = MultiBasisFunctionsState!(C, N, %s, %s);", filterSpec.hasOption!"usePower", buildArgs()));
}


/*
final class GeneralPowerState(C, size_t N, size_t P, size_t Mf = 0, size_t Mp = 0, ubyte filterSpec = 0)
if(P > 1)
{
  private
  {
    enum bool withIQImbalance = (filterSpec & FilterOptions.withConjugate) != 0;
    enum bool withDCBias = (filterSpec & FilterOptions.withDCBias) != 0;
    enum bool useConjugate = (filterSpec & FilterOptions.useConjugate) != 0;
    enum bool usePower = (filterSpec & FilterOptions.usePower) != 0;

    static assert(!withDCBias);
  }


    alias F = typeof(C.re);

    C[Mf+1] ymemory;
    C[Mf+1+Mp] xmemory;
    F[Mf+1+Mp] pmemory;
    C[(Mf+1+Mp) * (withIQImbalance ? 2 : 1)][N] state;
    C[(Mf+1+Mp) * (withIQImbalance ? 2 : 1)][N] weight;
    size_t cnt;

  static if(usePower)
  {
    F[(Mf+1+Mp) * (withIQImbalance ? 2 : 1)] power;
  }


    this(F initP)
    {
        foreach(ref e; xmemory) e = complexZero!C;
        foreach(ref e; ymemory) e = complexZero!C;
        foreach(ref e; pmemory) e = 0;
        foreach(ref e; state) foreach(ref ee; e) ee = complexZero!C;
        foreach(ref e; weight) foreach(ref ee; e) ee = complexZero!C;
        static if(usePower) foreach(ref e; power) e = initP;
    }


    private enum string updateMethod =
    q{
        xmemory.shiftBack;
        ymemory.shiftBack;
        pmemory.shiftBack;
        state.shiftBack;

        xmemory[0] = c;
        immutable cabs = c.abs();

        pmemory[0] = cabs^^(P-1);

        immutable xk = xmemory[Mf];

        foreach(m; 0 .. Mp+Mf+1){
            auto xx = xk * pmemory[m];

            state[0][m] = xx;

          static if(withIQImbalance)
            state[0][(Mf+Mp+1) + m] = xx.conj;
        }

      static if(usePower)
        foreach(i; 0 .. power.length){
            auto e = state[0][i];
            power[i] = e.re^^2 + e.im^^2;
        }
    };


    private enum string errorMethod =
    q{
        ymemory[0] = c;
        c = ymemory[Mf];
        foreach(i; 0 .. state.length) foreach(j; 0 .. state[0].length){
            c -= state[i][j] * weight[i][j];
        }
    };


    void update(C c) pure nothrow @safe @nogc
    {
        mixin(updateMethod);
    }


    C error(C c) pure nothrow @safe @nogc
    {
        mixin(errorMethod);
        return c;
    }


    void apply(in C[] txs, in C[] rxs, C[] dst) pure nothrow @safe @nogc
    {
        foreach(ic, tx; txs){
            C c = tx;
            static if(useConjugate) c = c.conj;
            mixin(updateMethod);

            c = rxs[ic];
            mixin(errorMethod);

            dst[ic] = c;
        }
    }
}*/


final class BiasState(C)
{
    alias F = typeof(C.re);

    C[1][1] state;
    C[1][1] weight;
    F[1] power;


    this()
    {
        power[0] = 1;
        state[0][0] = 1 + complexZero!C;
        weight[0][0] = complexZero!C;
    }


    void update(C c) {}
    C error(C c)
    {
        return c - state[0][0] * weight[0][0];
    }

    void apply(in C[] tx, in C[] rx, C[] dst)
    {
        auto bias = state[0][0] * weight[0][0];

        foreach(i, e; tx)
            dst[i] = rx[i] - bias;
    }
}


auto inputTransformer(alias f, State, T...)(State state, T args)
{
    return new InputTransformer!(State, f, T)(state, args);
}


final class InputTransformer(S, alias f, T...)
{
    this(S state, T args)
    {
        sp = state;
        this.args = args;
    }


    void update(C)(C c)
    {
        sp.update(f(c, args));
    }


    C error(C)(C c)
    {
        return sp.error(c);
    }


    void apply(C)(in C[] tx, in C[] rx, C[] dst)
    {
        foreach(ic, C c; tx)
        {
            this.update(f(c, args));
            dst[ic] = this.error(rx[ic]);
        }
    }


    S sp;
    T args;

    alias sp this;
}


final class OneTapMultiFIRFilterState(C, size_t P)
{
    C[P][1] state;
    C[P][1] weight;
    typeof(C.re)[P] power;


    this()
    {
        foreach(ref e; state[0]) e = complexZero!C;
        foreach(ref e; weight[0]) e = complexZero!C;
        foreach(ref e; power) e = 1;
    }


    void update(in ref C[P] x)
    {
        state[0] = x;
        foreach(i, ref e; power)power[0] = x[i].sqAbs;
    }


    C error(C y)
    {
        foreach(i; 0 .. P)
            y -= state[0][i] * weight[0][i];

        return y;
    }
}
