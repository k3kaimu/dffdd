module dffdd.filter.lms;

import std.complex;
import std.math;
import std.range : lockstep;

import mir.ndslice;


final class NLMSAdapter(State)
{
    import std.algorithm;

    alias C = State.StateElementType;
    alias F = typeof(C.init.re);


    this(State state, real mu)
    {
        _mu = mu;
    }


    void adapt(ref State state, C error) /*pure nothrow @safe @nogc*/
    {
        real power = 0;
        foreach(j; 0 .. state.numOfFIR)
        foreach(i; 0 .. state.numOfTaps)
            power += state.state[i, j].sqAbs;

        foreach(j; 0 .. state.numOfFIR){
            foreach(i; 0 .. state.numOfTaps)
                state.weight[i, j] += _mu * error * state.state[i, j].conj / power;
        }
    }


  private:
    immutable real _mu;
}


auto makeNLMSAdapter(State)(State state, real mu)
{
    return new NLMSAdapter!State(state, mu);
}

