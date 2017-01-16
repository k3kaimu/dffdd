module dffdd.filter.filter;

import std.algorithm;
import std.complex;
import std.math;
import std.typetuple;
import std.range;
import std.traits;
import std.experimental.ndslice;
import std.numeric;

import std.stdio;

import carbon.math;
import carbon.stream;
import carbon.traits;

import dffdd.filter.state;
import dffdd.filter.traits;
import dffdd.utils.fft;
import dffdd.filter.primitives;


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
    auto state1 = MultiFIRState!(Complex!float)(1, 1),
         state2 = MultiFIRState!(Complex!float)(1, 1);

    // FIRフィルタとLMSアルゴリズムでフィルタを構成
    auto filter = parallelFilter(
        state1, makeLMSAdapter(state1, 0.1),
        state2, makeLMSAdapter(state2, 0.1)
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

/+
/**
C: Complex Type
numOfFIRState: the number of FIR filters
func[0]: function(C[]) -> C[numOfFIRState][], distortion function
func[1]: function(MultiFIRState) -> Adapter
*/
final class GeneralParallelHammersteinFilter(C, size_t numOfFIRState, func...)
if(func.length == 2)
{
    enum bool isFunc0Type = is(func[0]);    // isType!(func[0]) と等価

  static if(isFunc0Type)
  {
    this(func[0] distorter, size_t numOfTaps)
    {
        _f = distorter;
        _state = MultiFIRState!C(numOfFIRState, numOfTaps);
        _adapter = func[1](_state);
    }
  }
  else
  {
    this(size_t numOfTaps)
    {
        _state = MultiFIRState!C(numOfFIRState, numOfTaps);
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
        auto distorted = func[0](tx);

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
    MultiFIRState!C _state;
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
    auto filter = generalParallelHammersteinFilter!(Cpx, N, distortionFunc, s => makeLMSAdapter(s, 0.1))(10);

    // こんな感じでも生成できる
    auto filter2 = generalParallelHammersteinFilter!(Cpx, N, s => makeLMSAdapter(s, 0.1))
                    (
                        delegate (Cpx[] tx) => tx.map!(a => cast(Cpx[N])[a, a*(a.abs^^2)]),
                        10
                    );
}
+/



final class SimpleTimeDomainParallelHammersteinFilter(C, Dist, alias genAdaptor)
if(isBlockConverter!(Dist, C, C[]))
{
    this(Dist dist, size_t numOfFIR, size_t numOfTaps)
    in{
        assert(numOfFIR == dist.outputDim);
    }
    body {
        _distorter = dist;
        _state = MultiFIRState!(Complex!float)(numOfFIR, numOfTaps);
        _adaptor = genAdaptor(_state);
        _buffer = new C[][](this.inputBlockLength, dist.outputDim);
    }


    size_t inputBlockLength() @property
    {
        return _distorter.inputBlockLength;
    }


    void apply(Flag!"learning" doLearning)(in C[] input, in C[] desired, C[] errors)
    in{
        assert(input.length == desired.length);
        assert(input.length == errors.length);
        assert(input.length % this.inputBlockLength == 0);
    }
    body{
        immutable size_t blk = this.inputBlockLength;

        foreach(i; 0 .. input.length / blk){
            auto ips = input[i*blk .. (i+1) * blk];
            auto dss = desired[i*blk .. (i+1) * blk];
            auto ers = errors[i*blk .. (i+1) * blk];

            _distorter(ips, _buffer);
            foreach(j, e; _buffer){
                _state.update(e);
                ers[j] = dss[j] - _state.output;

              static if(doLearning)
                _adaptor.adapt(_state, ers[j]);
            }
        }
    }


  static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
  {
    private void learningFromTX(R)(R txs)
    if(isInputRange!R && is(Unqual!(ElementType!R) : C))
    {
        _distorter.learn(txs);
    }
  }


    void preLearning(R1, R2, R3)(R1 digitalTx, R2 paDirects, R3 exampleSI)
    {
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            learningFromTX(digitalTx);
    }


  private:
    Dist _distorter;
    MultiFIRState!(Complex!float) _state;
    typeof(genAdaptor(_state)) _adaptor;
    C[][] _buffer;
}

final class SimpleTimeDomainParallelHammersteinFilterWithDCM(C, Dist, alias genAdaptor)
if(isBlockConverter!(Dist, C, C[]))
{
    this(Dist dist, size_t numOfFIR, size_t numOfTaps)
    in{
        assert(numOfFIR == dist.outputDim);
    }
    body {
        _distorter = dist;
        _state = MultiFIRState!(Complex!float)(numOfFIR, numOfTaps);

        foreach(i; 0 .. numOfFIR)
            _adaptors ~= genAdaptor(_state.subFIRState[i]);

        _buffer = new C[][](this.inputBlockLength, this.outputDim);
    }


    size_t inputBlockLength() @property
    {
        return _distorter.inputBlockLength;
    }


    void apply(Flag!"learning" doLearning)(in C[] input, in C[] desired, C[] errors)
    in{
        assert(input.length == desired.length);
        assert(input.length == errors.length);
        assert(input.length % this.inputBlockLength == 0);
    }
    body{
        immutable size_t blk = this.inputBlockLength;

        foreach(i; 0 .. input.length / blk){
            auto ips = input[i*blk .. (i+1) * blk];
            auto dss = desired[i*blk .. (i+1) * blk];
            auto ers = errors[i*blk .. (i+1) * blk];

            _distorter(ips, _buffer);
            foreach(j, e; _buffer){
                _state.update(e);
                ers[j] = dss[j] - _state.output;

              static if(doLearning)
              {
                foreach(k, ref e; _adaptors){
                    auto subState = _state.subFIRState(k);
                    e.adapt(subState, ers[j]);
                }
              }
            }
        }
    }


    Dist distorter() @property { return _distorter; }


  static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
  {
    private void learningFromTX(R)(R txs)
    if(isInputRange!R && is(Unqual!(ElementType!R) : C))
    {
        _distorter.learn(txs);
    }
  }


    void preLearning(R1, R2, R3)(R1 digitalTx, R2 paDirects, R3 exampleSI)
    {
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            learningFromTX(digitalTx);
    }



  private:
    Dist _distorter;
    MultiFIRState!(Complex!float) _state;
    typeof(genAdaptor(_state))[] _adaptors;
    C[][] _buffer;
}


final class SimpleTimeDomainCascadeHammersteinFilter(C, Dist, alias genAdaptor)
if(isBlockConverter!(Dist, C, C[]))
{
    this(Dist dist, size_t numOfFIR, size_t numOfTaps)
    in{
        assert(numOfFIR == dist.outputDim);
    }
    body {
        _distorter = dist;
        _state = MultiFIRState!(Complex!float)(numOfFIR, numOfTaps);

        foreach(i; 0 .. numOfFIR)
            _adaptors ~= genAdaptor(i, _state.subFIRState(i));

        _buffer = new C[][](this.inputBlockLength, dist.outputDim);
    }


    size_t inputBlockLength() @property
    {
        return _distorter.inputBlockLength;
    }


    void apply(Flag!"learning" doLearning)(in C[] input, in C[] desired, C[] errors)
    in{
        assert(input.length == desired.length);
        assert(input.length == errors.length);
        assert(input.length % this.inputBlockLength == 0);
    }
    body{
        immutable size_t blk = this.inputBlockLength;

        foreach(i; 0 .. input.length / blk){
            auto ips = input[i*blk .. (i+1) * blk];
            auto dss = desired[i*blk .. (i+1) * blk];
            auto ers = errors[i*blk .. (i+1) * blk];

            _distorter(ips, _buffer);
            foreach(j, e; _buffer){
                _state.update(e);

                ers[j] = dss[j];
                foreach(k, ref adaptor; _adaptors){
                    auto subState = _state.subFIRState(k);
                    ers[j] -= subState.output;

                  static if(doLearning){
                    adaptor.adapt(subState, ers[j]);
                  }
                }
            }
        }
    }


    Dist distorter() @property { return _distorter; }


  static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
  {
    private void learningFromTX(R)(R txs)
    if(isInputRange!R && is(Unqual!(ElementType!R) : C))
    {
        _distorter.learn(txs);
    }
  }


    void preLearning(R1, R2, R3)(R1 digitalTx, R2 paDirects, R3 exampleSI)
    {
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            learningFromTX(digitalTx);
    }

  private:
    Dist _distorter;
    MultiFIRState!(Complex!float) _state;
    typeof(genAdaptor(0, _state))[] _adaptors;
    C[][] _buffer;
}


final class SimpleFrequencyDomainParallelHammersteinFilter(C, Dist, SpecConv, alias genAdaptor)
if(isBlockConverter!(Dist, C, C[]))
{
    enum bool doesDistHaveLearn = is(typeof((Dist dist, C[] r){ dist.learn(r); }));

    import dffdd.filter.regenerator : OverlapSaveRegenerator;

    this(Dist dist, SpecConv specConv, in bool[] subcarrierMap, size_t nFFT, size_t nCP, size_t nOS)
    {
        immutable dim = dist.outputDim;

        _fftw = makeFFTWObject!Complex(nFFT * nOS);
        _distorter = dist;
        _specConv = specConv;
        _invSpecConv = specConv;
        auto sliceState = new C[nFFT * nOS * dim].sliced(nFFT * nOS, dim);
        sliceState[] = complexZero!C;

        _regenerator = new OverlapSaveRegenerator!C(sliceState);

        foreach(i; 0 .. nFFT * nOS){
            auto state = MultiFIRState!(Complex!float)(dim, 1);
            // state.weight = sliceState[i .. (i+1), 0 .. $];
            _states ~= state;
            _adaptors ~= genAdaptor(i, subcarrierMap[i], state);
        }

        _scMap = subcarrierMap.dup;
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;
        _distortedBuffer = new C[][](this.inputBlockLength, dim);
        _fftedBuffer = new C[][](dim, _nFFT * _nOS);
    }


    size_t inputBlockLength() @property
    {
        immutable a = _distorter.inputBlockLength,
                  b = (_nFFT + _nCP) * _nOS;

        return a * (b / gcd(a, b));     // LCM
    }


    void apply(Flag!"learning" doLearning)(in C[] input, in C[] desired, C[] errors)
    in{
        assert(input.length == desired.length);
        assert(input.length == errors.length);
        assert(input.length % this.inputBlockLength == 0);
        foreach(e; input){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
        foreach(e; desired){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    out{
        foreach(e; errors){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    body{
        immutable size_t blk = this.inputBlockLength;
        immutable size_t dim = _distorter.outputDim;

        if(doLearning)
        {
            foreach(i; 0 .. input.length / blk){
                auto ips = input[i*blk .. (i+1) * blk];
                auto dss = desired[i*blk .. (i+1) * blk];
                auto ers = errors[i*blk .. (i+1) * blk];

                _distorter(ips, _distortedBuffer);
                
                immutable symLen = _nOS * (_nFFT + _nCP),
                          cpLen = _nOS * _nCP;

                foreach(j; 0 .. blk / symLen){
                    auto bsym = _distortedBuffer[j * symLen + cpLen .. (j + 1) * symLen];
                    auto dsym = dss[j * symLen + cpLen .. (j + 1) * symLen];

                    // 非線形送信信号を周波数領域に変換
                    foreach(p; 0 .. dim){
                        foreach(k; 0 .. _nFFT * _nOS)
                            _fftw.inputs!float[k] = bsym[k][p];

                        _fftw.fft!float();
                        _fftedBuffer[p][] = _fftw.outputs!float[];
                    }

                    // 受信信号の周波数成分
                    _fftw.inputs!float[] = dsym[];
                    _fftw.fft!float();
                    auto rxFFTed = _fftw.outputs!float();

                    // 全周波数で，適応フィルタにかける
                    C[] vX = new C[dim];
                    C[] vZ = new C[dim];
                    foreach(f; 0 .. _nFFT * _nOS){
                        foreach(p, ref e; vX)
                            e = _fftedBuffer[p][f];

                        // _specConvで送信信号のスペクトルに細工を加える
                      static if(!is(SpecConv == typeof(null)))
                        _specConv.converter(f)(vX, vZ);
                      else
                        vZ[] = vX[];

                        _states[f].update(vZ);
                        C error = rxFFTed[f] - _states[f].output;
                        _adaptors[f].adapt(_states[f], error);
                    }
                }
            }

            // 除去用の周波数応答へ逆変換する
            C[] vH = new C[dim];
            C[] vG = new C[dim];
            foreach(f; 0 .. _nFFT * _nOS){
                vH.sliced(dim)[] = _states[f].weight[0];

              static if(!is(SpecConv == typeof(null)))
                _invSpecConv.converter(f)(vH, vG);
              else
                vG[] = vH[];

                _regenerator.frequencyResponse[f, 0 .. $] = vG.sliced(dim);
            }
        }

        // 除去
        foreach(i; 0 .. input.length / blk){
            auto ips = input[i*blk .. (i+1) * blk];
            auto dss = desired[i*blk .. (i+1) * blk];
            auto ers = errors[i*blk .. (i+1) * blk];

            _distorter(ips, _distortedBuffer);
            foreach(es; _distortedBuffer) foreach(e; es){
                assert(!isNaN(e.re));
                assert(!isNaN(e.im));
            }

            immutable symLen = _nOS * (_nFFT + _nCP),
                      cpLen = _nOS * _nCP;

            _regenerator.apply(_distortedBuffer, ers);
            foreach(e; ers){
                assert(!isNaN(e.re));
                assert(!isNaN(e.im));
            }

            foreach(j; 0 .. blk)
                ers[j] = dss[j] - ers[j];
        }
    }


    Dist distorter() @property { return _distorter; }


    private void learningFromTX(R)(R txs)
    if(isInputRange!R && is(Unqual!(ElementType!R) : C))
    {
        immutable size_t dim = _distorter.outputDim;

      static if(doesDistHaveLearn)
      {
        _distorter.learn(txs.save);
      }

      static if(!is(SpecConv == typeof(null)))
      {
        // TODO
        C[][] xs;
        foreach(e; txs){
            xs.length += 1;
            _distorter(e, xs[$-1]);
        }

        immutable symLen = _nOS * (_nFFT + _nCP),
                  cpLen = _nOS * _nCP;

        C[][][] trs;
        foreach(i; 0 .. xs.length / symLen){
            C[][] spec = new C[][](_nFFT * _nOS, dim);
            auto sym = xs[i * symLen + cpLen .. (i + 1) * symLen];
            foreach(p; 0 .. dim)
            {
                auto ips = _fftw.inputs!float;
                auto ops = _fftw.outputs!float;
                foreach(j; 0 .. _nFFT * _nOS)
                    ips[j] = sym[j][p];
                
                _fftw.fft!float();
                foreach(j; 0 .. _nFFT * _nOS)
                    spec[j][p] = ops[j];
            }

            trs ~= spec;
        }

        _specConv.learn(trs, _scMap);
        import std.stdio;
        writeln(_specConv.converter(1).convertMatrix);
        _invSpecConv = _specConv.inverter;
        writeln(_invSpecConv.converter(1).convertMatrix);
      }
    }


    /**

    */
    void preLearning(R1, R2, R3)(R1 digitalTx, R2 paDirects, R3 exampleSI)
    {
        learningFromTX(digitalTx);
    }


  private:
    FFTWObject!Complex _fftw;
    Dist _distorter;
    SpecConv _specConv, _invSpecConv;
    OverlapSaveRegenerator!C _regenerator;
    immutable(bool[]) _scMap;
    immutable size_t _nFFT, _nCP, _nOS;
    MultiFIRState!(Complex!float)[] _states;
    typeof(genAdaptor(0, true, _states[0]))[] _adaptors;
    C[][] _distortedBuffer;
    C[][] _fftedBuffer;
}


final class SimpleFrequencyDomainCascadeHammersteinFilter(C, Dist, SpecConv, alias genAdaptor)
if(isBlockConverter!(Dist, C, C[]))
{
    enum bool doesDistHaveLearn = is(typeof((Dist dist, C[] r){ dist.learn(r); }));

    import dffdd.filter.regenerator : OverlapSaveRegenerator;

    this(Dist dist, SpecConv specConv, in bool[] subcarrierMap, size_t nFFT, size_t nCP, size_t nOS)
    {
        immutable dim = dist.outputDim;

        _fftw = makeFFTWObject!Complex(nFFT * nOS);
        _distorter = dist;
        _specConv = specConv;
        _invSpecConv = specConv;
        auto sliceState = new C[nFFT * nOS * dim].sliced(nFFT * nOS, dim);
        sliceState[] = complexZero!C;

        _regenerator = new OverlapSaveRegenerator!C(sliceState);

        foreach(i; 0 .. nFFT * nOS){
            MultiFIRState!(Complex!float)[] sts;
            typeof(_adaptors[0]) ads;
            foreach(p; 0 .. dim){
                auto state = MultiFIRState!(Complex!float)(1, 1);
                // state.weight = sliceState[i .. (i+1), 0 .. $];
                sts ~= state;
                ads ~= genAdaptor(i, subcarrierMap[i], p, state);
            }
            _states ~= sts;
            _adaptors ~= ads;
        }

        _scMap = subcarrierMap.dup;
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;
        _distortedBuffer = new C[][](this.inputBlockLength, dim);
        _fftedBuffer = new C[][](dim, _nFFT * _nOS);
    }


    size_t inputBlockLength() @property
    {
        immutable a = _distorter.inputBlockLength,
                  b = (_nFFT + _nCP) * _nOS;

        return a * (b / gcd(a, b));     // LCM
    }


    void apply(Flag!"learning" doLearning)(in C[] input, in C[] desired, C[] errors)
    in{
        assert(input.length == desired.length);
        assert(input.length == errors.length);
        assert(input.length % this.inputBlockLength == 0);
        foreach(e; input){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
        foreach(e; desired){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    out{
        foreach(e; errors){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    body{
        immutable size_t blk = this.inputBlockLength;
        immutable size_t dim = _distorter.outputDim;

        if(doLearning)
        {
            foreach(i; 0 .. input.length / blk){
                auto ips = input[i*blk .. (i+1) * blk];
                auto dss = desired[i*blk .. (i+1) * blk];
                auto ers = errors[i*blk .. (i+1) * blk];

                _distorter(ips, _distortedBuffer);
                
                immutable symLen = _nOS * (_nFFT + _nCP),
                          cpLen = _nOS * _nCP;

                foreach(j; 0 .. blk / symLen){
                    auto bsym = _distortedBuffer[j * symLen + cpLen .. (j + 1) * symLen];
                    auto dsym = dss[j * symLen + cpLen .. (j + 1) * symLen];

                    // 非線形送信信号を周波数領域に変換
                    foreach(p; 0 .. dim){
                        foreach(k; 0 .. _nFFT * _nOS)
                            _fftw.inputs!float[k] = bsym[k][p];

                        _fftw.fft!float();
                        _fftedBuffer[p][] = _fftw.outputs!float[];
                    }

                    // 受信信号の周波数成分
                    _fftw.inputs!float[] = dsym[];
                    _fftw.fft!float();
                    auto rxFFTed = _fftw.outputs!float();

                    // 全周波数で，適応フィルタにかける
                    C[] vX = new C[dim];
                    C[] vZ = new C[dim];
                    foreach(f; 0 .. _nFFT * _nOS){
                        foreach(p, ref e; vX)
                            e = _fftedBuffer[p][f];

                        // _specConvで送信信号のスペクトルに細工を加える
                      static if(!is(SpecConv == typeof(null)))
                        _specConv.converter(f)(vX, vZ);
                      else
                        vZ[] = vX[];

                        C error = rxFFTed[f];
                        foreach(p; 0 .. dim){
                            _states[f][p].update(vZ[p]);
                            error -= _states[f][p].output;
                            _adaptors[f][p].adapt(_states[f][p], error);
                        }
                    }
                }
            }

            // 除去用の周波数応答へ逆変換する
            C[] vH = new C[dim];
            C[] vG = new C[dim];
            foreach(f; 0 .. _nFFT * _nOS){
                foreach(p; 0 .. dim)
                    vH[p] = _states[f][p].weight[0, 0];

              static if(!is(SpecConv == typeof(null)))
              {
                _invSpecConv.converter(f)(vH, vG);
                foreach(p; 0 .. dim)
                    _regenerator.frequencyResponse[f, p] = vG[p];
              }
              else
              {
                foreach(p; 0 .. dim)
                    _regenerator.frequencyResponse[f, p] = vH[p];
              }
            }
        }

        // 除去
        foreach(i; 0 .. input.length / blk){
            auto ips = input[i*blk .. (i+1) * blk];
            auto dss = desired[i*blk .. (i+1) * blk];
            auto ers = errors[i*blk .. (i+1) * blk];

            _distorter(ips, _distortedBuffer);
            foreach(es; _distortedBuffer) foreach(e; es){
                assert(!isNaN(e.re));
                assert(!isNaN(e.im));
            }

            _regenerator.apply(_distortedBuffer, ers);
            foreach(e; ers){
                assert(!isNaN(e.re));
                assert(!isNaN(e.im));
            }

            foreach(j; 0 .. blk)
                ers[j] = dss[j] - ers[j];
        }
    }


    Dist distorter() @property { return _distorter; }


  //   static if(is(typeof((SpecConv specConv, C[][][] txs){ specConv.learn(txs); })))
  //   {
    void learningFromTX(R)(R txs)
    if(isInputRange!R && is(Unqual!(ElementType!R) : C))
    {
        immutable size_t dim = _distorter.outputDim;

      static if(doesDistHaveLearn)
      {
        _distorter.learn(txs);
      }

      static if(!is(SpecConv == typeof(null)))
      {
        // TODO
        C[][] xs;
        foreach(e; txs){
            xs.length += 1;
            _distorter(e, xs[$-1]);
        }

        immutable symLen = _nOS * (_nFFT + _nCP),
                  cpLen = _nOS * _nCP;

        C[][][] trs;
        foreach(i; 0 .. xs.length / symLen){
            C[][] spec = new C[][](_nFFT * _nOS, dim);
            auto sym = xs[i * symLen + cpLen .. (i + 1) * symLen];
            foreach(p; 0 .. dim)
            {
                auto ips = _fftw.inputs!float;
                auto ops = _fftw.outputs!float;
                foreach(j; 0 .. _nFFT * _nOS)
                    ips[j] = sym[j][p];
                
                _fftw.fft!float();
                foreach(j; 0 .. _nFFT * _nOS)
                    spec[j][p] = ops[j];
            }

            trs ~= spec;
        }

        _specConv.learn(trs, _scMap);
        import std.stdio;
        writeln(_specConv.converter(1).convertMatrix);
        _invSpecConv = _specConv.inverter;
        writeln(_invSpecConv.converter(1).convertMatrix);
        // _invSpecConv = _specConv;

        // import std.stdio;
        // _distorter.learn(txs);
        // writeln(_distorter.converter.convertMatrix);
      }
    }
  //   }


    void preLearning(R1, R2, R3)(R1 digitalTx, R2 paDirects, R3 exampleSI)
    {
        learningFromTX(digitalTx.map!(x => cast(C)x));
    }


  private:
    FFTWObject!Complex _fftw;
    Dist _distorter;
    SpecConv _specConv, _invSpecConv;
    OverlapSaveRegenerator!C _regenerator;
    immutable(bool[]) _scMap;
    immutable size_t _nFFT, _nCP, _nOS;
    MultiFIRState!(Complex!float)[][] _states;
    typeof(genAdaptor(0, true, 0, _states[0][0]))[][] _adaptors;
    C[][] _distortedBuffer;
    C[][] _fftedBuffer;
}


final class SimpleFrequencyDomainCascadeHammersteinFilterType2(C, Dist, alias genAdaptor)
{
    enum bool doesDistHaveLearn = is(typeof((Dist dist, C[] r){ dist.learn(r); }));


    this(Dist dist, in bool[] subcarrierMap, size_t nFFT, size_t nCP, size_t nOS)
    {
        _distorter = dist;

        foreach(p; aliasSeqOf!(iota(0, Dist.outputDim))){
            _filters[p] = new typeof(_filters[p])(new LinearDist(), null, subcarrierMap, nFFT, nCP, nOS);

            foreach(i; 0 .. nFFT * nOS){
                auto s = _filters[p]._states[i];
                _filters[p]._adaptors[i] = genAdaptor(i, subcarrierMap[i], p, s);
            }
        }

        _inputBuf = new C[][](Dist.outputDim, this.inputBlockLength);
        _errorBuf = new C[](this.inputBlockLength);
    }


    size_t inputBlockLength() @property
    {
        return _filters[0].inputBlockLength;
    }


    void apply(Flag!"learning" doLearning)(in C[] input, in C[] desired, C[] errors)
    in{
        assert(input.length % this.inputBlockLength == 0);
        assert(input.length == desired.length);
        assert(input.length == errors.length);
    }
    body{
        immutable blen = this.inputBlockLength;

        foreach(i; 0 .. input.length / blen){

            auto iblock = input[i*blen .. (i+1)*blen];
            auto dblock = desired[i*blen .. (i+1)*blen];
            auto eblock = errors[i*blen .. (i+1)*blen];

            // 直交化処理
            foreach(j; 0 .. blen)
            {
                C[Dist.outputDim] outs_;
                auto tmpBuf = outs_[];
                _distorter(iblock[j], tmpBuf);
                foreach(p; 0 .. Dist.outputDim)
                    _inputBuf[p][j] = tmpBuf[p];
            }

            _errorBuf[] = dblock[];

            foreach(p; aliasSeqOf!(iota(0, Dist.outputDim))){
                _filters[p].apply!(doLearning)(_inputBuf[p], _errorBuf, eblock);
                _errorBuf[] = eblock[];
            }
        }
    }


    Dist distorter() @property { return _distorter; }


    private void learningFromTX(R)(R txs)
    if(isInputRange!R && is(Unqual!(ElementType!R) : C))
    {
        immutable size_t dim = _distorter.outputDim;

      static if(doesDistHaveLearn)
      {
        _distorter.learn(txs.save);
      }
    }


    void preLearning(R1, R2, R3)(R1 digitalTx, R2 paDirects, R3 exampleSI)
    {
        learningFromTX(digitalTx.map!(x => cast(C)x));
    }


  private:
    Dist _distorter;
    FilterTypes!(Dist.outputDim) _filters;
    C[][] _inputBuf;
    C[] _errorBuf;

  static
  {
    // template genAdaptorImpl(size_t p)
    // {
    //     auto genAdaptorImpl(State)(size_t fIdx, bool isUsedSC, State state)
    //     {
    //         return ;
    //     }
    // }

    alias LinearDist = Distorter!(C, x => x);

    template FilterTypes(size_t P)
    {
      static if(P == 0)
        alias FilterTypes = AliasSeq!();
      else
        alias FilterTypes = AliasSeq!(FilterTypes!(P-1),
                                      SimpleFrequencyDomainParallelHammersteinFilter!(C, LinearDist, typeof(null), (f, b, s) => typeof(genAdaptor(f, b, P-1, s)).init));
    }
  }
}
