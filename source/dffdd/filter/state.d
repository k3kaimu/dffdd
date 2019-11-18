module dffdd.filter.state;

import std.complex;
import mir.ndslice;

import carbon.math;


struct MultiFIRState(C, SliceKind kind = Contiguous)
if(kind == Contiguous || kind == Canonical)
{
    alias StateElementType = C;
    alias R = typeof(C.init.re);

    Slice!(C*, 2, kind) state, weight;

    this(size_t numOfFIR, size_t nTaps)
    {
        static if(kind == Contiguous)
        {
            this(numOfFIR, nTaps, 
                new C[nTaps * numOfFIR].sliced(nTaps, numOfFIR),
                new C[nTaps * numOfFIR].sliced(nTaps, numOfFIR));
        }
        else
        {
            this(numOfFIR, nTaps, 
                new C[nTaps * numOfFIR].sliced(nTaps, numOfFIR).canonical,
                new C[nTaps * numOfFIR].sliced(nTaps, numOfFIR).canonical);
        }

        this.state[] = C(0);
        this.weight[] = C(0);
    }


    this(size_t numOfFIR, size_t nTaps, Slice!(C*, 2, kind) state, Slice!(C*, 2, kind) weight)
    in{
        assert(state.length!0 == nTaps && weight.length!0 == nTaps);
    }
    body
    {
        this.state = state;
        this.weight = weight;
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
    }


    C output() @property
    {
        import dffdd.utils.linalg : dot;
        C dst = C(0);

        static if(kind == Contiguous)
        {
            dst += dot(state.flattened, weight.flattened);
        }
        else
        {
            foreach(i; 0 .. this.numOfTaps)
                dst += dot(state[i], weight[i]);
        }

        return dst;
    }


    MultiFIRState!(C, Canonical) subFIRState(size_t i)
    {
        return MultiFIRState!(C, Canonical)(1, this.numOfTaps, state[0 .. $, i .. i+1], weight[0 .. $, i .. i+1]);
    }


    MultiFIRState!(C, Canonical) subFIRStates(size_t i, size_t j)
    {
        return MultiFIRState!(C, Canonical)(1, this.numOfTaps, state[0 .. $, i .. j], weight[0 .. $, i .. j]);
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
