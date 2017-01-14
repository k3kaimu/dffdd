module dffdd.filter.state;

import std.complex;
import std.experimental.ndslice;

import carbon.math;




//final class FIRState(C, bool usePower = false)
//{
//    alias StateElementType = C;
//    //enum bool hasMultiStateFIRFilter = true;
//    //enum size_t numOfInputs = 1;
//    //enum size_t numOfOutputs = 1;
//    enum size_t numOfFIRState = 1;
//    //enum size_t numOfFIRTaps = NumOfTaps;

//    Slice!(1, C*) state, weight;

//  static if(usePower)
//    typeof(C.init.re) power;


//    this(size_t nTaps)
//    {
//        this(nTaps, new C[nTaps].sliced(nTaps),
//                    new C[nTaps].sliced(nTaps));
//    }


//    this(size_t nTaps, Slice!(1, C*) state, Slice!(1, C*) weight)
//    {
//        this.state = state;
//        this.weight = weight;

//        this.state[] = complexZero!C;
//        this.weight[] = complexZero!C;

//      static if(usePower)
//        this.power = 1;
//    }


//    size_t numOfTaps() const @property { return state.length!0; }


//    void update(C x)
//    {
//        foreach_reverse(i; 1 .. this.numOfTaps)
//            state[i] = state[i-1];

//        state[0] = x;
//        static if(usePower) power = x.re ^^2 + x.im ^^ 2;
//    }


//    C error(C y)
//    {
//        foreach(i; 0 .. this.numOfTaps)
//            y -= state[i] * weight[i];

//        return y;
//    }
//}


//auto inputTransformer(alias f, State, T...)(State state, T args)
//{
//    return new InputTransformer!(State, f, T)(state, args);
//}



//final class InputTransformer(S, alias f, T...)
//{
//    alias StateElementType = S.StateElementType;
//    //enum bool hasMultiFIRState = S.hasMultiFIRState;
//    //enum size_t numOfFIRState = S.numOfFIRState;
//    //enum size_t numOfFIRTaps = S;

//    alias numOfFIRTaps = sp.numOfFIRState;
//    alias numOfFIRState = sp.numOfFIRState;


//    this(S state, T args)
//    {
//        sp = state;
//        this.args = args;
//    }


//    void update(C)(C c)
//    {
//        sp.update(f(c, args));
//    }


//    C error(C)(C c)
//    {
//        return sp.error(c);
//    }


//    void apply(C)(in C[] tx, in C[] rx, C[] dst)
//    {
//        foreach(ic, C c; tx)
//        {
//            this.update(f(c, args));
//            dst[ic] = this.error(rx[ic]);
//        }
//    }


//    S sp;
//    T args;

//    alias sp this;
//}



//final class MultiFIRState(C, size_t NumOfFIRState, bool usePower = false)
//{
//    alias StateElementType = C;
//    enum size_t numOfFIRState = NumOfFIRState;

//    Slice!(2, C*) state, weight;
//    static if(usePower) Slice!(1, typeof(C.init.re)*) power;

//    this(size_t nTaps)
//    {
//        this(nTaps, new C[nTaps * NumOfFIRState].sliced(nTaps, NumOfFIRState),
//                    new C[nTaps * NumOfFIRState].sliced(nTaps, NumOfFIRState));
//    }


//    this(size_t nTaps, Slice!(2, C*) state, Slice!(2, C*) weight)
//    {
//        this.state = state;
//        this.weight = weight;

//        this.state[] = complexZero!C;
//        this.weight[] = complexZero!C;

//      static if(usePower)
//      {
//        this.power = new typeof(C.init.re)[NumOfFIRState].sliced(NumOfFIRState);
//        this.power[] = 1;
//      }
//    }



//    size_t numOfTaps() const @property { return state.length!0; }


//    void update(C[NumOfFIRState] x...)
//    {
//        foreach_reverse(i; 1 .. this.numOfTaps)
//            state[i][] = state[i-1][];

//        state[0][] = x[];
//        static if(usePower)
//            foreach(i; 0 .. NumOfFIRState)
//                power[i] = x[i].re ^^2 + x[i].im ^^ 2;
//    }


//    C error(C y)
//    {
//        foreach(i; 0 .. this.numOfTaps)
//            foreach(j; 0 .. this.numOfFIRState)
//                y -= state[i, j] * weight[i, j];

//        return y;
//    }
//}


//final class ParallelHammersteinState(C, size_t NumOfBasisFuncs, bool usePower = false)
//{
//    MultiFIRState!(C, NumOfBasisFuncs, usePower) multiFIRState;
//    alias multiFIRState this;


//    this(size_t nTaps, C delegate(C)[NumOfBasisFuncs] funcs...)
//    {
//        multiFIRState = new MultiFIRState!(C, NumOfBasisFuncs, true)(nTaps);
//        _funcs = funcs;
//    }


//    void update(C x)
//    {
//        C[NumOfBasisFuncs] vs;
//        foreach(i; 0 .. NumOfBasisFuncs)
//            vs[i] = _funcs[i](x);

//        multiFIRState.update(vs);
//    }


//  private:
//    C delegate(C)[NumOfBasisFuncs] _funcs;
//}


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

struct MultiFIRState(C)
{
    alias StateElementType = C;
    alias R = typeof(C.init.re);

    Slice!(2, C*) state, weight;
    Slice!(1, R*) power;

    this(size_t numOfFIR, size_t nTaps)
    {
        this(numOfFIR, nTaps, 
            new C[nTaps * numOfFIR].sliced(nTaps, numOfFIR),
            new C[nTaps * numOfFIR].sliced(nTaps, numOfFIR),
            new R[numOfFIR].sliced(numOfFIR));

        this.state[] = complexZero!C;
        this.weight[] = complexZero!C;
        this.power[] = 1;
    }


    this(size_t numOfFIR, size_t nTaps, Slice!(2, C*) state, Slice!(2, C*) weight, Slice!(1, R*) power)
    in{
        assert(state.length!0 == nTaps && weight.length!0 == nTaps);
        assert(state.length!1 == numOfFIR && weight.length!1 == numOfFIR && power.length == numOfFIR);
    }
    body
    {
        this.state = state;
        this.weight = weight;
        this.power = power;
    }


    size_t numOfTaps() const @property { return state.length!0; }
    size_t numOfFIR() const @property { return state.length!1; }


    void update(C[] x...)
    in{
        assert(x.length == numOfFIR);
    }
    body{
        foreach_reverse(i; 1 .. state.length!0)
            state[i][] = state[i-1][];

        state[0][] = x[];
        foreach(i; 0 .. state.length!1)
            power[i] = x[i].re ^^2 + x[i].im ^^ 2;
    }


    C output() @property
    {
        C dst = complexZero!C;

        foreach(i; 0 .. this.numOfTaps)
            foreach(j; 0 .. this.numOfFIR)
                dst += state[i, j] * weight[i, j];

        return dst;
    }


    MultiFIRState!C subFIRState(size_t i)
    {
        return MultiFIRState!C(1, this.numOfTaps, state[0 .. $, i .. i+1], weight[0 .. $, i .. i+1], power[i .. i+1]);
    }


    MultiFIRState!C subFIRStates(size_t i, size_t j)
    {
        return MultiFIRState!C(1, this.numOfTaps, state[0 .. $, i .. j], weight[0 .. $, i .. j], power[i .. j]);
    }
}

unittest
{
    import std.math;

    auto fir = MultiFIRState!(Complex!float)(1, 1);

    assert(fir.state.shape == [1, 1]);
    assert(fir.state[0, 0] == 0);
    assert(fir.weight[0, 0] == 0);
    assert(fir.power[0] == 1);

    fir.update(Complex!float(1, 1));
    assert(fir.state[0, 0] == Complex!float(1, 1));
    assert(fir.weight[0, 0] == 0);
    assert(fir.power[0].approxEqual(2));
    assert(fir.output == 0);

    fir.weight[0, 0] = 1;
    assert(fir.output == Complex!float(1, 1));
    fir.weight[0, 0] = 2;
    assert(fir.output == Complex!float(2, 2));

    fir = MultiFIRState!(Complex!float)(2, 3);

    assert(fir.state.shape == [3, 2]);
    foreach(i; 0 .. 3) foreach(j; 0 .. 2) {
        assert(fir.state[i, j] == 0);
        assert(fir.weight[i, j] == 0);
        assert(fir.power[j] == 1);
    }

    fir.weight[0, 0] = 1;
    fir.weight[0, 1] = 2;
    fir.weight[1, 0] = 3;
    fir.weight[1, 1] = 4;

    fir.update(Complex!float(1, 1), Complex!float(1, 0));
    assert(fir.state[0, 0] == Complex!float(1, 1));
    assert(fir.state[0, 1] == Complex!float(1, 0));
    assert(fir.state[1, 0] == 0);
    assert(fir.state[1, 1] == 0);
    assert(fir.output == 1 * Complex!float(1, 1) + 2 * Complex!float(1, 0));

    fir.update(Complex!float(0, 1), Complex!float(-1, -1));
    assert(fir.state[0, 0] == Complex!float(0, 1));
    assert(fir.state[0, 1] == Complex!float(-1, -1));
    assert(fir.state[1, 0] == Complex!float(1, 1));
    assert(fir.state[1, 1] == Complex!float(1, 0));

    assert(fir.output == 1 * Complex!float(0, 1) + 2 * Complex!float(-1, -1)
                       + 3 * Complex!float(1, 1) + 4 * Complex!float(1, 0));
}
