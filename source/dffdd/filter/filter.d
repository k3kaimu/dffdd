module dffdd.filter.filter;

import std.algorithm;
import std.complex;
import std.math;
import std.typetuple;
import std.range;
import std.traits;
import std.experimental.ndslice;
import std.numeric;
import std.typecons;

import std.stdio;

import carbon.math;
import carbon.stream;
import carbon.traits;

import dffdd.blockdiagram.model;

import dffdd.filter.state;
import dffdd.filter.traits;
import dffdd.utils.fft;
import dffdd.filter.primitives;
import dffdd.utils.unit;
import dffdd.dsp.statistics;


template hasMemberAsSignal(S, names...)
if(names.length == 0 || is(typeof(names[0]) : string))
{
  static if(names.length == 0)
  {
      enum bool hasMemberAsSignal = true;
  }
  else
  {
      enum bool hasMemberAsSignal = is(typeof(mixin("S.init." ~ names[0]))) && isInputRange!(typeof(mixin("S.init." ~ names[0]))) && hasMemberAsSignal!(S, names[1 .. $]);
  }
}


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


//   static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
//   {
//     private void learningFromTX(R)(R txs)
//     if(isInputRange!R && is(Unqual!(ElementType!R) : C))
//     {
//         _distorter.learn(txs);
//     }
//   }


//     void preLearning(R1, R2, R3)(R1 digitalTx, R2 paDirects, R3 exampleSI)
//     {
//         static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
//             learningFromTX(digitalTx);
//     }


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator)
    if(isModelParameterSet!M)
    {
        // メンバーを持っているかどうかチェック
        static assert(hasMemberAsSignal!(typeof(signalGenerator(model)), "txBaseband"));

        // learningFromTX(signalGenerator(model).txBaseband);
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            _distorter.learn(signalGenerator(model).txBaseband);
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
            _adapters ~= genAdaptor(_state.subFIRState[i]);

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
                foreach(k, ref e; _adapters){
                    auto subState = _state.subFIRState(k);
                    e.adapt(subState, ers[j]);
                }
              }
            }
        }
    }


    Dist distorter() @property { return _distorter; }


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator)
    if(isModelParameterSet!M)
    {
        // メンバーを持っているかどうかチェック
        static assert(hasMemberAsSignal!(typeof(signalGenerator(model)), "txBaseband"));

        // learningFromTX(signalGenerator(model).txBaseband);
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            _distorter.learn(signalGenerator(model).txBaseband);
    }



  private:
    Dist _distorter;
    MultiFIRState!(Complex!float) _state;
    typeof(genAdaptor(_state))[] _adapters;
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
            _adapters ~= genAdaptor(i, _state.subFIRState(i));

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
                foreach(k, ref adaptor; _adapters){
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


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator)
    if(isModelParameterSet!M)
    {
        // メンバーを持っているかどうかチェック
        static assert(hasMemberAsSignal!(typeof(signalGenerator(model)), "txBaseband"));

        // learningFromTX(signalGenerator(model).txBaseband);
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            _distorter.learn(signalGenerator(model).txBaseband);
    }


  private:
    Dist _distorter;
    MultiFIRState!(Complex!float) _state;
    typeof(genAdaptor(0, _state))[] _adapters;
    C[][] _buffer;
}


final class SimpleFrequencyDomainParallelHammersteinFilter(C, Dist, Adapter, bool doProfiling = false)
if(isBlockConverter!(Dist, C, C[]))
{
    enum bool doesDistHaveLearn = is(typeof((Dist dist, C[] r){ dist.learn(r); }));

    import dffdd.filter.regenerator : OverlapSaveRegenerator;

    this(Dist dist, Adapter delegate(size_t i, bool b, MultiFIRState!C) genAdapter, in bool[] subcarrierMap, size_t nFFT, size_t nCP, size_t nOS, real sampFreq, Gain limitIRR = 20.dB, Gain gainMargin = 10.dB)
    {
        immutable dim = dist.outputDim;

        _fftw = makeFFTWObject!Complex(nFFT * nOS);
        _distorter = dist;
        _genAdapter = genAdapter;

        auto sliceState = new C[nFFT * nOS * dim].sliced(nFFT * nOS, dim);
        sliceState[] = complexZero!C;

        _regenerator = new OverlapSaveRegenerator!C(sliceState);

        foreach(i; 0 .. nFFT * nOS){
            auto state = MultiFIRState!(Complex!float)(dim, 1);
            // state.weight = sliceState[i .. (i+1), 0 .. $];
            _states ~= state;
            _adapters ~= genAdapter(i, subcarrierMap[i], state);
        }

        _scMap = subcarrierMap.dup;
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;
        _sampFreq = sampFreq;
        _distortedBuffer = new C[][](this.inputBlockLength, dim);
        _fftedBuffer = new C[][](dim, _nFFT * _nOS);
        _limitIRR = limitIRR;
        _gainMargin = gainMargin;
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
            if(_importance !is null && _selectedBasisFuncs is null){
                // 
                _selectedBasisFuncs = new bool[][](_nFFT * _nOS, Dist.outputDim);

                auto txs = input[_nCP * _nOS .. (_nFFT + _nCP) * _nOS].dup;
                auto rxs = desired[_nCP * _nOS .. (_nFFT + _nCP) * _nOS].dup;

                auto fftw = globalBankOf!(makeFFTWObject)[_nFFT * _nOS];
                fftw.inputs!float[] = txs[];
                fftw.fft!float();
                txs[] = fftw.outputs!float[];

                fftw.inputs!float[] = rxs[];
                fftw.fft!float();
                rxs[] = fftw.outputs!float[];

                foreach(f; 0 .. _nFFT * _nOS){
                    float h = (rxs[f] * txs[f].conj / (txs[f].sqAbs + _noiseFloor)).sqAbs;
                    float hinv = (rxs[$-1-f] * txs[$-1-f].conj / (txs[$-1-f].sqAbs + _noiseFloor)).sqAbs;

                    writefln("%s, %s, %s", f, h, hinv);

                    // この条件が満たされたとき，全基底関数を使用する
                    if(hinv > h * _limitIRR.gain^^2) {
                        // writefln("%s, %s, %s", f, h, hinv);
                        foreach(p; 0 .. Dist.outputDim)
                            _selectedBasisFuncs[f][p] = true;
                        
                        continue;
                    }

                    size_t cnt;
                    foreach(p; 0 .. Dist.outputDim) if(!_selectedBasisFuncs[f][p]) {
                        // 重要度とチャネルの積が，推定された干渉電力になる
                        float ip = _importance[f][p] * rxs[f].sqAbs;
                        // if(f == 0 || f == 256){
                        //     writefln("%s, %s, %s, %s, %s", p, ip, _importance[f][p], rxs[f].sqAbs, _noiseFloor);
                        // }

                        // 干渉電力にマージンを設けて，ノイズ電力と比較する
                        if(ip * _gainMargin.gain^^2 > _noiseFloor){
                            _selectedBasisFuncs[f][p] = true;
                            _selectedBasisFuncs[f][Dist.indexOfConjugated(p)] = true;
                            cnt += 2;
                        }
                    }

                    if(cnt == 0){
                        // writeln(_importance[f]);
                        // writeln(h);
                        // writeln(hinv);
                        // writeln(rxs[f].sqAbs);
                        // writeln(_noiseFloor);
                        // writeln(_importance[f][0] * rxs[f].sqAbs);
                        // assert(cnt != 0);
                        _selectedBasisFuncs[f][0] = true;
                        cnt += 1;
                    }

                    // _states, _adaptersを更新する
                    _states[f] = MultiFIRState!(Complex!float)(cnt, 1);
                    _adapters[f] = _genAdapter(f, _scMap[f], _states[f]);
                }

                // foreach(f; 0 .. _nFFT * _nOS)
                //     writefln("%s, %(%s, %)", f, _selectedBasisFuncs[f].map!`a ? 1 : 0`.zip(iota(1, _selectedBasisFuncs[f].length+1)).map!`a[0]*a[1]`);
            }

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
                    foreach(f; 0 .. _nFFT * _nOS){
                        foreach(p, ref e; vX)
                            e = _fftedBuffer[p][f];

                        auto vX2 = vX;
                        if(_selectedBasisFuncs !is null) {
                            // vXの中身を，必要なものだけにする
                            size_t cnt = 0;
                            foreach(p; 0 .. dim)
                                if(_selectedBasisFuncs[f][p]){
                                    // 選ばれているなら
                                    vX[cnt] = vX[p];
                                    cnt += 1;
                                }
                            
                            vX2 = vX[0 .. cnt];
                        }

                        _states[f].update(vX2);
                        C error = rxFFTed[f] - _states[f].output;
                        _adapters[f].adapt(_states[f], error);
                    }
                }
            }

            // 学習した結果をRegeneratorに設定する
            foreach(f; 0 .. _nFFT * _nOS){
                if(_selectedBasisFuncs !is null) {
                    size_t cnt = 0;
                    foreach(p; 0 .. dim){
                        if(_selectedBasisFuncs[f][p]){
                            _regenerator.frequencyResponse[f, p] = _states[f].weight[0, cnt];
                            cnt += 1;
                        }else{
                            _regenerator.frequencyResponse[f, p] = cast(C)Complex!float(0, 0);
                        }
                    }
                }else{
                    _regenerator.frequencyResponse[f, 0 .. $] = _states[f].weight[0];
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


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator)
    // if(isModelParameterSet!M)
    {
        // メンバーを持っているかどうかチェック
        static assert(hasMemberAsSignal!(typeof(signalGenerator(model)), "txBaseband", "noise", "receivedSISWP"));

        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            _distorter.learn(signalGenerator(model).txBaseband);

        static if(doProfiling)
            profiling(model, signalGenerator);
    }


    void profiling(M, Signals)(M model, Signals delegate(M) genSignal)
    {
        import dffdd.filter.ls : LSAdapter, makeLSAdapter;

        // 同軸線路を使用するように設定する
        model.useCoaxialCableAsChannel();

        // 最初にノイズ電力の計測
        // _noiseFloor = signals.noise.save.map!`a.re^^2 + a.im^^2`.take(1024*1024).sum() / (1024 * 1024);
        {
            auto signals = genSignal(model);

            auto n = signals.noise.save;
            auto fftw = globalBankOf!(makeFFTWObject)[_nFFT * _nOS];
            float[] buf = new float[_nFFT * _nOS];
        
            foreach(ref e; buf)
                e = 0;

            // auto buf = new Complex!float[_nOS * _nFFT];
            foreach(i; 0 .. 1024){
                foreach(j; 0 .. _nOS * _nFFT){
                    fftw.inputs!float[j] = cast(Complex!float)n.front;
                    n.popFront();
                }

                fftw.fft!float();
                foreach(j; 0 .. _nOS * _nFFT)
                    buf[j] += fftw.outputs!float[j].sqAbs;
            }

            foreach(j; 0 .. _nOS * _nFFT)
                buf[j] /= 1024;

            _noiseFloor = buf.sum() / (_nFFT * _nOS);
        }

        // LNAダイナミックレンジ最大まで使用する
        model.INR = model.lna.DR;

        // {
        //     auto model4IRR = model;
        //     model.ofdm.numOfSubcarrier = 1;
        //     model.noImpairementsOnTX();
        // }

        auto signals = genSignal(model);

        // 次のようなフィルタを学習してみる
        // + 適応アルゴリズム : 最小二乗法
        // + 使用シンボル数：1024シンボル
        auto testfilter = new SimpleFrequencyDomainParallelHammersteinFilter!(C, Dist, LSAdapter!(MultiFIRState!C))
                        (   _distorter,
                            (size_t i, bool b, MultiFIRState!C s) => makeLSAdapter(s, 1024),
                            _scMap, _nFFT, _nCP, _nOS, _sampFreq);

        immutable spb = this.inputBlockLength / ((_nFFT + _nCP) * _nOS);    // 1ブロックあたりのシンボル数

        C[] testTX = new C[this.inputBlockLength],
            testRX = new C[this.inputBlockLength],
            testER = new C[this.inputBlockLength];

        // 1024シンボル使用して学習する
        foreach(i; 0 .. 1024 / spb){
            signals.fillBuffer!(["txBasebandSWP", "receivedSISWP"])(testTX, testRX);
            testfilter.apply!(Yes.learning)(testTX, testRX, testER);
        }

        // 各周波数で，学習完了後の重みを得る
        C[][] weight = new C[][](_nFFT * _nOS, _distorter.outputDim);
        foreach(f; 0 .. _nFFT * _nOS)
            foreach(p; 0 .. _distorter.outputDim)
                weight[f][p] = testfilter._states[f].weight[0][p];

        real[][] powerOfSI = new real[][](_nFFT * _nOS, _distorter.outputDim);
        real[][] firstPowerOfSI = new real[][](_nFFT * _nOS, _distorter.outputDim);
        real[] powerOfRCV = new real[](_nFFT * _nOS);
        real[] firstPowerOfRCV = new real[](_nFFT * _nOS);
        foreach(f; 0 .. _nFFT * _nOS){
            powerOfRCV[f] = 0;
            foreach(p; 0 .. _distorter.outputDim)
                powerOfSI[f][p] = 0;
        }

        testTX = testTX[0 .. (_nFFT + _nCP) * _nOS];
        testRX = testRX[0 .. (_nFFT + _nCP) * _nOS];
        auto disted = new C[][]((_nFFT + _nCP) * _nOS, _distorter.outputDim);

        // 重要度を計算する
        foreach(i; 0 .. 1024){
            signals.fillBuffer!(["txBasebandSWP", "receivedSISWP"])(testTX, testRX);
            _distorter(testTX, disted);

            // 非線形送信信号を周波数領域に変換
            foreach(p; 0 .. _distorter.outputDim){
                foreach(k; 0 .. _nFFT * _nOS)
                    _fftw.inputs!float[k] = disted[_nCP * _nOS + k][p];

                _fftw.fft!float();

                foreach(k; 0 .. _nFFT * _nOS){
                    if(i == 0) firstPowerOfSI[k][p] = _fftw.outputs!float[k].sqAbs;
                    disted[k][p] = _fftw.outputs!float[k] * weight[k][p];
                }
            }

            // 受信信号の周波数成分
            _fftw.inputs!float[] = testRX[_nCP * _nOS .. $];
            _fftw.fft!float();
            auto rxFFTed = _fftw.outputs!float();

            // 電力に足し合わせる
            foreach(f; 0 .. _nFFT * _nOS){
                foreach(p; 0 .. _distorter.outputDim){
                    powerOfSI[f][p] += disted[f][p].sqAbs;
                }

                powerOfRCV[f] += rxFFTed[f].sqAbs;
            }

            if(i == 0)
                firstPowerOfRCV = powerOfRCV.dup;
        }

        _importance = powerOfSI;
        foreach(f; 0 .. _nFFT * _nOS) foreach(p; 0 .. _distorter.outputDim){
            _importance[f][p] /= powerOfRCV[f];
        }

        // auto flipped = _importance[$/2 .. $].chain(_importance[0 .. $/2]);
        // auto prcvFlped = powerOfRCV[$/2 .. $].chain(powerOfRCV[0 .. $/2]);
        // auto fpsFlped = firstPowerOfSI[$/2 .. $].chain(firstPowerOfSI[0 .. $/2]);
        // auto fprFlped = firstPowerOfRCV[$/2 .. $].chain(firstPowerOfRCV[0 .. $/2]);
        // foreach(f; 0 .. _nFFT * _nOS){
        //     writefln("%s, %(%s, %), %s, %s, %s, %s, %s", f, flipped[f].map!(a => 10*log10(a)), 10*log10(prcvFlped[f]), 10*log10(prcvFlped[_nFFT * _nOS -1 - f]),
        //         10*log10(fprFlped[f] / fpsFlped[f][0]),
        //         10*log10(fprFlped[$-1-f] / fpsFlped[$-1-f][0]),
        //         10*log10(fprFlped[$-1-f] / fpsFlped[$-1-f][0]) - 10*log10(fprFlped[f] / fpsFlped[f][0]) > 20 ? 0 : 100,
        //     );
        // }
    }


  private:
    FFTWObject!Complex _fftw;
    Dist _distorter;
    OverlapSaveRegenerator!C _regenerator;
    Adapter delegate(size_t, bool, MultiFIRState!C) _genAdapter;
    immutable(bool[]) _scMap;
    immutable size_t _nFFT, _nCP, _nOS;
    immutable real _sampFreq;
    immutable Gain _limitIRR, _gainMargin;
    MultiFIRState!C[] _states;
    Adapter[] _adapters;
    C[][] _distortedBuffer;
    C[][] _fftedBuffer;
    real[][] _importance;
    bool[][] _selectedBasisFuncs;
    float _noiseFloor = 0;
}


final class SimpleFrequencyDomainDCMHammersteinFilterType2(C, Dist, alias genAdaptor, Flag!"isParallel" isParallel = No.isParallel)
if(isBlockConverter!(Dist, C, C[]))
{
    enum bool doesDistHaveLearn = is(typeof((Dist dist, C[] r){ dist.learn(r); }));

    import dffdd.filter.regenerator : OverlapSaveRegenerator;

    this(Dist dist, in bool[] subcarrierMap, size_t nFFT, size_t nCP, size_t nOS)
    {
        immutable dim = dist.outputDim;

        _fftw = makeFFTWObject!Complex(nFFT * nOS);
        _distorter = dist;
        auto sliceState = new C[nFFT * nOS * dim].sliced(nFFT * nOS, dim);
        sliceState[] = complexZero!C;

        _regenerator = new OverlapSaveRegenerator!C(sliceState);

        foreach(i; 0 .. nFFT * nOS){
            MultiFIRState!(Complex!float)[] sts;
            typeof(_adapters[0]) ads;
            foreach(p; 0 .. dim){
                auto state = MultiFIRState!(Complex!float)(1, 1);
                // state.weight = sliceState[i .. (i+1), 0 .. $];
                sts ~= state;
                ads ~= genAdaptor(i, subcarrierMap[i], p, state);
            }
            _states ~= sts;
            _adapters ~= ads;
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
                    foreach(f; 0 .. _nFFT * _nOS){
                        foreach(p, ref e; vX)
                            e = _fftedBuffer[p][f];

                        C error = rxFFTed[f];
                        foreach(p; 0 .. dim){
                            _states[f][p].update(vX[p]);
                            error -= _states[f][p].output;

                          static if(!isParallel)
                            _adapters[f][p].adapt(_states[f][p], error);
                        }

                      static if(isParallel)
                      {
                        foreach(p; 0 .. dim)
                            _adapters[f][p].adapt(_states[f][p], error);
                      }
                    }
                }
            }

            // 学習した結果をregeneratorに伝達する
            C[] vH = new C[dim];
            foreach(f; 0 .. _nFFT * _nOS){
                foreach(p; 0 .. dim)
                    vH[p] = _states[f][p].weight[0, 0];

                foreach(p; 0 .. dim)
                    _regenerator.frequencyResponse[f, p] = vH[p];
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


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator)
    if(isModelParameterSet!M)
    {
        // メンバーを持っているかどうかチェック
        static assert(hasMemberAsSignal!(typeof(signalGenerator(model)), "txBaseband"));

        // learningFromTX(signalGenerator(model).txBaseband);
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            _distorter.learn(signalGenerator(model).txBaseband);
    }


  private:
    FFTWObject!Complex _fftw;
    Dist _distorter;
    OverlapSaveRegenerator!C _regenerator;
    immutable(bool[]) _scMap;
    immutable size_t _nFFT, _nCP, _nOS;
    MultiFIRState!(Complex!float)[][] _states;
    typeof(genAdaptor(0, true, 0, _states[0][0]))[][] _adapters;
    C[][] _distortedBuffer;
    C[][] _fftedBuffer;
}


final class SimpleFrequencyDomainDCMHammersteinFilterType1(C, Dist, alias genAdaptor, Flag!"isParallel" isParallel = No.isParallel)
{
    enum bool doesDistHaveLearn = is(typeof((Dist dist, C[] r){ dist.learn(r); }));

    this(Dist dist, in bool[] subcarrierMap, size_t nFFT, size_t nCP, size_t nOS)
    {
        _distorter = dist;

        foreach(p; aliasSeqOf!(iota(0, Dist.outputDim))){
            _filters[p] = new typeof(_filters[p])(new LinearDist(), subcarrierMap, nFFT, nCP, nOS);

            foreach(i; 0 .. nFFT * nOS){
                auto s = _filters[p]._states[i];
                _filters[p]._adapters[i] = genAdaptor(i, subcarrierMap[i], p, s);
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


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator)
    if(isModelParameterSet!M)
    {
        // メンバーを持っているかどうかチェック
        static assert(hasMemberAsSignal!(typeof(signalGenerator(model)), "txBaseband"));

        // learningFromTX(signalGenerator(model).txBaseband);
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            _distorter.learn(signalGenerator(model).txBaseband);
    }



  private:
    Dist _distorter;
    FilterTypes!(Dist.outputDim) _filters;
    C[][] _inputBuf;
    C[] _errorBuf;

  static
  {
    alias LinearDist = Distorter!(C, x => x);

    template FilterTypes(size_t P)
    {
      static if(P == 0)
        alias FilterTypes = AliasSeq!();
      else
        alias FilterTypes = AliasSeq!(FilterTypes!(P-1),
                                      SimpleFrequencyDomainParallelHammersteinFilter!(C, LinearDist, (f, b, s) => typeof(genAdaptor(f, b, P-1, s)).init));
    }
  }
}
