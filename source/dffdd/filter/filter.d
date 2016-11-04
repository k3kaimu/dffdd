module dffdd.filter.polynomial;

import std.algorithm;
import std.complex;
import std.math;
import std.typetuple;

import std.stdio;

import carbon.stream;
import carbon.traits;

import dffdd.filter.traits;


void learning(F, C)(ref F filter, in C[] tx, in C[] rx, C[] outputBuf)
{
    filter.apply!true(tx, tx, outputBuf);
}


void cancelling(F, C)(ref F filter, in C[] tx, in C[] rx, C[] outputBuf)
{
    filter.apply!false(tx, rx, outputBuf);
}


final class PolynomialFilter(State, Adapter)
{
    alias C = typeof(State.init.state[0][0]);
    size_t cnt;

    this(State state, Adapter adapter)
    {
        _state = state;
        _adapter = adapter;
    }


    void apply(bool bLearning = true, C1 : C)(in C1[] tx, in C1[] rx, C1[] outputBuf)
    in{
        assert(tx.length == rx.length
            && tx.length == outputBuf.length);
    }
    body{
        foreach(i; 0 .. tx.length){
            _state.update(tx[i]);

            C error = _state.error(rx[i]);
            outputBuf[i] = error;

          static if(bLearning)
            _adapter.adapt(_state, error);
        }
        cnt += tx.length;
    }


    State state() @property { return _state; }


  private:
    State _state;
    Adapter _adapter;
}


auto polynomialFilter(State, Adapter)(State state, Adapter adapter)
{
    return new PolynomialFilter!(State, Adapter)(state, adapter);
}


/**
複数のフィルタが、直列に連結された状態のフィルタを構築します。
*/
final class SerialFilter(StateORAdapter...)
if(StateORAdapter.length % 2 == 0)
{
    template GetStride2(size_t n)
    {
      static if(n >= StateORAdapter.length)
        alias GetStride2 = TypeTuple!();
      else
        alias GetStride2 = TypeTuple!(StateORAdapter[n], GetStride2!(n+2));
    }

    alias States = GetStride2!(0);
    alias Adapters = GetStride2!(1);

    static assert(States.length == StateORAdapter.length / 2);
    static assert(Adapters.length == StateORAdapter.length / 2);

    States states;
    Adapters adapters;

    //alias C = typeof(states[0].state[0][0]);
    alias C = States[0].StateElementType;


    this(StateORAdapter stateAndAdapter)
    {
        foreach(i, T; StateORAdapter)
        {
          static if(i % 2 == 0)
            states[i/2] = stateAndAdapter[i];
          else
            adapters[i/2] = stateAndAdapter[i];
        }
    }


    void apply(bool bLearning, C1 : C)(in C1[] tx, in C1[] rx, C1[] outputBuf)
    in{
        assert(tx.length == rx.length
            && tx.length == outputBuf.length);
    }
    body{
        foreach(i; 0 .. tx.length){
            C error = rx[i];
            C txv = tx[i];
            foreach(j, T; States){
                states[j].update(txv);
                error = states[j].error(error);

              static if(bLearning)
                adapters[j].adapt(states[j], error);
            }
            outputBuf[i] = error;
        }
    }
}


auto serialFilter(StateORAdapter...)(StateORAdapter stateAndAdapter)
if(StateORAdapter.length % 2 == 0)
{
    return new SerialFilter!StateORAdapter(stateAndAdapter);
}

/**
複数のフィルタが、並列接続された状態のフィルタを構築します。
*/
final class ParallelFilter(StateORAdapter...)
if(StateORAdapter.length % 2 == 0)
{
    template GetStride2(size_t n)
    {
      static if(n >= StateORAdapter.length)
        alias GetStride2 = TypeTuple!();
      else
        alias GetStride2 = TypeTuple!(StateORAdapter[n], GetStride2!(n+2));
    }

    alias States = GetStride2!(0);
    alias Adapters = GetStride2!(1);

    States states;
    Adapters adaptors;

    alias C = States[0].StateElementType;


    this(StateORAdapter stateAndAdapter)
    {
        foreach(i, T; StateORAdapter)
        {
          static if(i % 2 == 0)
            states[i/2] = stateAndAdapter[i];
          else
            adaptors[i/2] = stateAndAdapter[i];
        }
    }


    void apply(bool bLearning = true, C1 : C)(in C1[] tx, in C1[] rx, C1[] outputBuf)
    in{
        assert(tx.length == rx.length
            && tx.length == outputBuf.length);
    }
    body{
        foreach(i; 0 .. tx.length){
            C error = rx[i];
            C txv = tx[i];
            foreach(j, T; States){
                states[j].update(txv);
                auto ev = states[j].error(rx[i]);
                error -= (rx[i] - ev);
            }

          static if(bLearning)
            foreach(j, T; States)
                adaptors[j].adapt(states[j], error);

            outputBuf[i] = error;
        }
    }
}


auto parallelFilter(StateORAdapter...)(StateORAdapter stateAndAdapter)
if(StateORAdapter.length % 2 == 0)
{
    return new ParallelFilter!StateORAdapter(stateAndAdapter);
}


/**
C: Complex Type
func: function(C[]) -> C[numOfFIRState][], distortion function
numOfFIRState: the number of FIR filters
genAdapter: function(MultiFIRState) -> Adapter
*/
final class GeneralParallelHammersteinFilter(C, alias func, size_t numOfFIRState, alias genAdapter)
{
    bool isFuncType = !is(typeof(func(null).front) : C[numOfFIRState]);

  static if(isFuncType)
  {
    this(func f, size_t numOfTaps)
    {
        _f = f;
        _state = new MultiFIRState!(C, numOfFIRState, true)(numOfTaps);
        _adapter = genAdapter(_state);
    }
  }
  else
  {
    this(size_t numOfTaps)
    {
        _state = new MultiFIRState!(C, numOfFIRState, true)(numOfTaps);
        _adapter = genAdapter(_state);
    }
  }


    void apply(bool bLearning = true, C1 : C)(in C1[] tx, in C1[] rx, C1[] outputBuf)
    in{
        assert(tx.length == rx.length
            && tx.length == outputBuf.length);
    }
    body{
      static if(isFuncType)
        auto distorted = _f(tx);
      else
        auto distorted = func(tx);

        foreach(i; 0 .. tx.length)
        {
            _state.update(distorted.front);
            C error = _state.error(rx[i]);
            outputBuf[i] = error;

          static if(bLearning)
            _adapter.adapt(_state, error);

            distorted.popFront();
        }
    }


  private:
    MultiFIRState!(C, numOfFIRState, true) _state;
    typeof(genAdapter(state)) _adapter;

  static if(isFuncType)
    func _f;
}
