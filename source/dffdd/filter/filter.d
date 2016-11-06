module dffdd.filter.polynomial;

import std.algorithm;
import std.complex;
import std.math;
import std.typetuple;
import std.range;
import std.traits;

import std.stdio;

import carbon.stream;
import carbon.traits;

import dffdd.filter.state;
import dffdd.filter.traits;


void learning(F, C)(ref F filter, in C[] tx, in C[] rx, C[] outputBuf)
{
    filter.apply!true(tx, tx, outputBuf);
}


void cancelling(F, C)(ref F filter, in C[] tx, in C[] rx, C[] outputBuf)
{
    filter.apply!false(tx, rx, outputBuf);
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


/**
複数のフィルタが並列に接続されたフィルタを構築します
*/
auto parallelFilter(StateORAdapter...)(StateORAdapter stateAndAdapter)
if(StateORAdapter.length % 2 == 0)
{
    return new ParallelFilter!StateORAdapter(stateAndAdapter);
}

///
unittest
{
    import dffdd.filter.state;
    import dffdd.filter.lms;

    // FIRフィルタ
    auto state1 = new FIRState!(Complex!float, true)(1),
         state2 = new FIRState!(Complex!float, true)(1);

    // FIRフィルタとLMSアルゴリズムでフィルタを構成
    auto filter = parallelFilter(
        state1, lmsAdapter(state1, 0.1),
        state2, lmsAdapter(state2, 0.1)
    );
}


/**
単一の状態と単一の適応アルゴリズムで構成されるフィルタを構築します．
parallelFilterの特殊版です．
*/
auto oneStateFilter(State, Adapter)(State state, Adapter adapter)
{
    return new ParallelFilter!(State, Adapter)(state, adapter);
}


/**
C: Complex Type
numOfFIRState: the number of FIR filters
func[0]: function(C[]) -> C[numOfFIRState][], distortion function
func[1]: function(MultiFIRState) -> Adapter
*/
final class GeneralParallelHammersteinFilter(C, size_t numOfFIRState, func...)
if(func.length == 2)
{
    enum bool isFunc0Type = isType!(func[0]);

  static if(isFunc0Type)
  {
    this(func[0] distorter, size_t numOfTaps)
    {
        _f = distorter;
        _state = new MultiFIRState!(C, numOfFIRState, true)(numOfTaps);
        _adapter = func[1](_state);
    }
  }
  else
  {
    this(size_t numOfTaps)
    {
        _state = new MultiFIRState!(C, numOfFIRState, true)(numOfTaps);
        _adapter = func[1](_state);
    }
  }


    void apply(bool bLearning = true, C1 : C)(in C1[] tx, in C1[] rx, C1[] outputBuf)
    in{
        assert(tx.length == rx.length
            && tx.length == outputBuf.length);
    }
    body{
      static if(isFunc0Type)
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
    typeof(func[1](_state)) _adapter;

  static if(isFunc0Type)
    func[0] _f;
}


/// ditto
auto generalParallelHammersteinFilter(C, size_t numOfFIRState, alias distorter, alias genAdapter)(size_t nTaps)
{
    return new GeneralParallelHammersteinFilter!(C, numOfFIRState, distorter, genAdapter)(nTaps);
}


/// ditto
auto generalParallelHammersteinFilter(C, size_t numOfFIRState, alias genAdapter, Func)(Func distorter, size_t nTaps)
{
    return new GeneralParallelHammersteinFilter!(C, numOfFIRState, Func, genAdapter)(distorter, nTaps);
}


///
unittest
{
    import dffdd.filter.state;
    import dffdd.filter.lms;
    import std.algorithm;
    import std.range;

    alias Cpx = Complex!float;
    enum size_t N = 2;  // 並列FIRフィルタの数

    // 入力を歪ませる関数
    static
    auto distortionFunc(Cpx[] tx)
    {
        return tx.map!(a => cast(Cpx[N])[a, a*(a.abs^^2)]);
    }

    // フィルタの生成
    auto filter = generalParallelHammersteinFilter!(Cpx, N, distortionFunc, s => lmsAdapter(s, 0.1))(10);

    // こんな感じでも生成できる
    auto filter2 = generalParallelHammersteinFilter!(Cpx, N, s => lmsAdapter(s, 0.1))
                    (
                        delegate (Cpx[] tx) => tx.map!(a => cast(Cpx[N])[a, a*(a.abs^^2)]),
                        10
                    );
}
