module dffdd.filter.state;

import std.complex;
import mir.ndslice;

import carbon.math;


struct MultiFIRState(C)
{
    alias StateElementType = C;
    alias R = typeof(C.init.re);

    Slice!(C*, 2, Universal) state, weight;
    // Slice!(Universal, [1], R*) power;

    this(size_t numOfFIR, size_t nTaps)
    {
        this(numOfFIR, nTaps, 
            new C[nTaps * numOfFIR].sliced(nTaps, numOfFIR).universal,
            new C[nTaps * numOfFIR].sliced(nTaps, numOfFIR).universal,
            /*new R[numOfFIR].sliced(numOfFIR).universal*/);

        this.state[] = complexZero!C;
        this.weight[] = complexZero!C;
        // this.power[] = 1;
    }


    this(size_t numOfFIR, size_t nTaps, Slice!(C*, 2, Universal) state, Slice!(C*, 2, Universal) weight/*, Slice!(Universal,[1], R*) power*/)
    in{
        assert(state.length!0 == nTaps && weight.length!0 == nTaps);
        // assert(state.length!1 == numOfFIR && weight.length!1 == numOfFIR && power.length == numOfFIR);
    }
    body
    {
        this.state = state;
        this.weight = weight;
        // this.power = power;
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
        // foreach(i; 0 .. state.length!1)
            // power[i] = x[i].re ^^2 + x[i].im ^^ 2;
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
        return MultiFIRState!C(1, this.numOfTaps, state[0 .. $, i .. i+1], weight[0 .. $, i .. i+1]/*, power[i .. i+1]*/);
    }


    MultiFIRState!C subFIRStates(size_t i, size_t j)
    {
        return MultiFIRState!C(1, this.numOfTaps, state[0 .. $, i .. j], weight[0 .. $, i .. j]/*, power[i .. j]*/);
    }
}

unittest
{
    import std.math;

    auto fir = MultiFIRState!(Complex!float)(1, 1);

    assert(fir.state.shape == [1, 1]);
    assert(fir.state[0, 0] == 0);
    assert(fir.weight[0, 0] == 0);
    // assert(fir.power[0] == 1);

    fir.update(Complex!float(1, 1));
    assert(fir.state[0, 0] == Complex!float(1, 1));
    assert(fir.weight[0, 0] == 0);
    // assert(fir.power[0].approxEqual(2));
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
        // assert(fir.power[j] == 1);
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
