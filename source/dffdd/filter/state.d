module dffdd.filter.state;

import std.complex;
import std.experimental.ndslice;

import carbon.math;




final class FIRState(C, bool usePower = false)
{
    alias StateElementType = C;
    //enum bool hasMultiStateFIRFilter = true;
    //enum size_t numOfInputs = 1;
    //enum size_t numOfOutputs = 1;
    enum size_t numOfFIRState = 1;
    //enum size_t numOfFIRTaps = NumOfTaps;

    Slice!(1, C*) state, weight;

  static if(usePower)
    typeof(C.init.re) power;


    this(size_t nTaps)
    {
        this(nTaps, new C[nTaps].sliced(nTaps),
                    new C[nTaps].sliced(nTaps));
    }


    this(size_t nTaps, Slice!(1, C*) state, Slice!(1, C*) weight)
    {
        this.state = state;
        this.weight = weight;

        this.state[] = complexZero!C;
        this.weight[] = complexZero!C;

      static if(usePower)
        this.power = 1;
    }


    size_t numOfTaps() const @property { return state.length!0; }


    void update(C x)
    {
        foreach_reverse(i; 1 .. this.numOfTaps)
            state[i] = state[i-1];

        state[0] = x;
        static if(usePower) power = x.re ^^2 + x.im ^^ 2;
    }


    C error(C y)
    {
        foreach(i; 0 .. this.numOfTaps)
            y -= state[i] * weight[i];

        return y;
    }
}


auto inputTransformer(alias f, State, T...)(State state, T args)
{
    return new InputTransformer!(State, f, T)(state, args);
}



final class InputTransformer(S, alias f, T...)
{
    alias StateElementType = S.StateElementType;
    //enum bool hasMultiFIRState = S.hasMultiFIRState;
    //enum size_t numOfFIRState = S.numOfFIRState;
    //enum size_t numOfFIRTaps = S;

    alias numOfFIRTaps = sp.numOfFIRState;
    alias numOfFIRState = sp.numOfFIRState;


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



final class MultiFIRState(C, size_t NumOfFIRState, bool usePower = false)
{
    alias StateElementType = C;
    enum size_t numOfFIRState = NumOfFIRState;

    Slice!(2, C*) state, weight;
    static if(usePower) Slice!(1, typeof(C.init.re)*) power;

    this(size_t nTaps)
    {
        this(nTaps, new C[nTaps * NumOfFIRState].sliced(nTaps, NumOfFIRState),
                    new C[nTaps * NumOfFIRState].sliced(nTaps, NumOfFIRState));
    }


    this(size_t nTaps, Slice!(2, C*) state, Slice!(2, C*) weight)
    {
        this.state = state;
        this.weight = weight;

        this.state[] = complexZero!C;
        this.weight[] = complexZero!C;

      static if(usePower)
      {
        this.power = new typeof(C.init.re)[NumOfFIRState].sliced(NumOfFIRState);
        this.power[] = 1;
      }
    }



    size_t numOfTaps() const @property { return state.length!0; }


    void update(C[NumOfFIRState] x...)
    {
        foreach_reverse(i; 1 .. this.numOfTaps)
            state[i][] = state[i-1][];

        state[0][] = x[];
        static if(usePower)
            foreach(i; 0 .. NumOfFIRState)
                power[i] = x[i].re ^^2 + x[i].im ^^ 2;
    }


    C error(C y)
    {
        foreach(i; 0 .. this.numOfTaps)
            foreach(j; 0 .. this.numOfFIRState)
                y -= state[i, j] * weight[i, j];

        return y;
    }
}


final class ParallelHammersteinState(C, size_t NumOfBasisFuncs, bool usePower = false)
{
    MultiFIRState!(C, NumOfBasisFuncs, usePower) multiFIRState;
    alias multiFIRState this;


    this(size_t nTaps, C delegate(C)[NumOfBasisFuncs] funcs...)
    {
        multiFIRState = new MultiFIRState!(C, NumOfBasisFuncs, true)(nTaps);
        _funcs = funcs;
    }


    void update(C x)
    {
        C[NumOfBasisFuncs] vs;
        foreach(i; 0 .. NumOfBasisFuncs)
            vs[i] = _funcs[i](x);

        multiFIRState.update(vs);
    }


  private:
    C delegate(C)[NumOfBasisFuncs] _funcs;
}



final class OneTapMultiFIRFilterState(C, size_t P)
{
    alias StateElementType = C;

    enum size_t numOfFIRState = P;
    Slice!(2, C*) state, weight;
    Slice!(1, typeof(C.init.re)*) power;


    this()
    {
        state = new C[P].sliced(1, P);
        weight = new C[P].sliced(1, P);
        power = new typeof(C.init.re)[P].sliced(P);

        foreach(ref e; state[0]) e = complexZero!C;
        foreach(ref e; weight[0]) e = complexZero!C;
        foreach(ref e; power) e = 1;
    }


    size_t numOfTaps() const @property { return state.length!0; }


    void update(in ref C[P] x)
    {
        state[0][] = x[];
        foreach(i; 0 .. P) power[i] = x[i].re^^2 + x[i].im^^2;
    }


    C error(C y)
    {
        foreach(i; 0 .. P)
            y -= state[0][i] * weight[0][i];

        return y;
    }
}

