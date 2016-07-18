module dffdd.filter.state;

import std.experimental.ndslice;

import carbon.math;




final class FIRFilter(C, bool usePower = false)
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
