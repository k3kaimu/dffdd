module dffdd.filter.filter;

import std.algorithm;
import std.complex;
import std.math;
import std.typetuple;
import std.range;
import std.traits;
import std.numeric;
import std.typecons;
import std.json;
import std.variant;

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
    do{
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
    do{
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
        state1, makeNLMSAdapter(state1, 0.1),
        state2, makeNLMSAdapter(state2, 0.1)
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


final class SimpleTimeDomainParallelHammersteinFilter(C, Dist, alias genAdaptor)
if(isBlockConverter!(Dist, C, C[]))
{
    this(Dist dist, size_t numOfTaps)
    {
        auto numOfFIR = dist.outputDim;
        _distorter = dist;
        _state = MultiFIRState!C(numOfFIR, numOfTaps);
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
    do{
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


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator)
    // if(isModelParameterSet!M)
    {
        // メンバーを持っているかどうかチェック
        static assert(hasMemberAsSignal!(typeof(signalGenerator(model)), "txBaseband"));

        // learningFromTX(signalGenerator(model).txBaseband);
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
        {
            auto sig = signalGenerator(model);
            auto buf = new C[](model.orthogonalizer.numOfTrainingSamples);
            sig.fillBuffer!(["txBaseband"])(buf);
            _distorter.learn(buf);
        }
    }


    ref MultiFIRState!C state() return @property 
    {
        return _state;
    }


  private:
    Dist _distorter;
    MultiFIRState!C _state;
    typeof(genAdaptor(_state)) _adaptor;
    C[][] _buffer;
}



final class SimpleTimeDomainMultiInputFilter(C, alias genAdaptor)
{
    this(size_t numInputs, size_t numOfTaps)
    {
        _state = MultiFIRState!C(numInputs, numOfTaps);
        _adaptor = genAdaptor(_state);
    }


    void apply(Flag!"learning" doLearning)(in C[][] inputs, in C[] desired, C[] errors)
    in{
        assert(inputs.length == desired.length);
        assert(inputs.length == errors.length);
        if(inputs.length > 0) assert(inputs[0].length == _state.numOfFIR);
    }
    do{
        foreach(i; 0 .. inputs.length) {
            _state.update(inputs[i]);
            errors[i] = desired[i] - _state.output;

            static if(doLearning)
                _adaptor.adapt(_state, errors[i]);
        }
    }


    ref MultiFIRState!C state() return @property 
    {
        return _state;
    }


  private:
    MultiFIRState!C _state;
    typeof(genAdaptor(_state)) _adaptor;
}



final class SimpleTimeDomainParallelHammersteinFilterWithDCM(C, Dist, alias genAdaptor)
if(isBlockConverter!(Dist, C, C[]))
{
    this(Dist dist, size_t numOfTaps)
    {
        auto numOfFIR = dist.outputDim;

        _distorter = dist;
        _state = MultiFIRState!C(numOfFIR, numOfTaps);

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
    do{
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
    {
        // メンバーを持っているかどうかチェック
        static assert(hasMemberAsSignal!(typeof(signalGenerator(model)), "txBaseband"));

        // learningFromTX(signalGenerator(model).txBaseband);
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            _distorter.learn(signalGenerator(model).txBaseband);
    }



  private:
    Dist _distorter;
    MultiFIRState!C _state;
    typeof(genAdaptor(_state))[] _adapters;
    C[][] _buffer;
}


final class SimpleTimeDomainCascadeHammersteinFilter(C, Dist, alias genAdaptor)
if(isBlockConverter!(Dist, C, C[]))
{
    this(Dist dist, size_t numOfTaps, size_t numIter = 2)
    {
        _distorter = dist;
        _state = MultiFIRState!C(dist.outputDim, numOfTaps);
        _numIter = numIter;

        foreach(i; 0 .. dist.outputDim)
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
    do{
        immutable size_t blk = this.inputBlockLength;

        if(doLearning) {
            _xsTr ~= input.dup;
            _ysTr ~= desired.dup;
            errors[] = desired[];
        } else {
            if(_lastState == State.train) {
                // 1ブランチずつ学習していく
                auto linear = new XApDistorter!(C, 1)();
                assert(linear.outputDim == 1);
                auto errorTr = new C[_ysTr.length];
                auto _distTrXs = new C[][](_xsTr.length, _distorter.outputDim);
                _distorter(_xsTr, _distTrXs);

                // _distTrXsを転置する
                auto distTrXs = new C[][](_distorter.outputDim, _xsTr.length);
                foreach(i; 0 .. _distorter.outputDim) foreach(j; 0 .. _xsTr.length)
                    distTrXs[i][j] = _distTrXs[j][i];

                auto lastRemained = _ysTr.dup;
                foreach(_; 0 .. this._numIter) {
                    foreach(i; 0 .. _distorter.outputDim) {
                        auto lincanc = new SimpleTimeDomainParallelHammersteinFilter!(C, typeof(linear), (s) => genAdaptor(i, s))(linear, this._state.numOfTaps);
                        lincanc.apply!(Yes.learning)(distTrXs[i], lastRemained, errorTr);
                        // 学習した結果で除去をして，次のブランチへ進む
                        lincanc.apply!(No.learning)(distTrXs[i], lastRemained, errorTr);
                        lastRemained[] = errorTr[];

                        // 学習した重みをコピーする
                        foreach(j; 0 .. this._state.numOfTaps) {
                            _state.weight[j, i] += lincanc.state.weight[j, 0];
                        }
                    }
                }

                // Stateを最新状態にしておく
                foreach(i; 0 .. _xsTr.length / blk){
                    auto ips = _xsTr[i*blk .. (i+1) * blk];
                    _distorter(ips, _buffer);
                    foreach(j, e; _buffer)
                        _state.update(e);
                }

                _xsTr = null;
                _ysTr = null;
            }

            // 除去する
            foreach(i; 0 .. input.length / blk){
                auto ips = input[i*blk .. (i+1) * blk];
                auto dss = desired[i*blk .. (i+1) * blk];
                auto ers = errors[i*blk .. (i+1) * blk];

                _distorter(ips, _buffer);
                foreach(j, e; _buffer){
                    _state.update(e);
                    ers[j] = dss[j] - _state.output;
                }
            }
        }

        _lastState = doLearning ? State.train : State.canc;
    }


    Dist distorter() @property { return _distorter; }


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator)
    {
        // メンバーを持っているかどうかチェック
        static assert(hasMemberAsSignal!(typeof(signalGenerator(model)), "txBaseband"));

        // learningFromTX(signalGenerator(model).txBaseband);
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            _distorter.learn(signalGenerator(model).txBaseband.map!(a => cast(C)a));
    }


  private:
    size_t _numIter;
    Dist _distorter;
    MultiFIRState!C _state;
    typeof(genAdaptor(0, _state.subFIRState(0)))[] _adapters;
    C[][] _buffer;
    C[] _xsTr, _ysTr;
    State _lastState = State.train;

    enum State { train, canc }
}


final class SimpleFrequencyDomainParallelHammersteinFilter(C, Dist, Adapter)
if(isBlockConverter!(Dist, C, C[]))
{
    alias R = typeof(C.init.re);
    enum bool doesDistHaveLearn = is(typeof((Dist dist, C[] r){ dist.learn(r); }));

    import dffdd.filter.regenerator : OverlapSaveRegenerator;

    this(Dist dist, Adapter delegate(size_t i, bool b, MultiFIRState!C) genAdapter, in bool[] subcarrierMap, size_t nFFT, size_t nCP, size_t nOS, real sampFreq, bool doSelect = false, size_t nEstH = 2, Gain limitIRR = 20.dB, Gain gainMargin = 6.dB)
    {
        immutable dim = dist.outputDim;

        _fftw = makeFFTWObject!Complex(nFFT * nOS);
        _distorter = dist;
        _genAdapter = genAdapter;

        auto sliceState = new C[nFFT * nOS * dim].sliced(nFFT * nOS, dim);
        sliceState[] = C(0);

        _regenerator = new OverlapSaveRegenerator!C(sliceState);

        foreach(i; 0 .. nFFT * nOS){
            auto state = MultiFIRState!C(dim, 1);
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
        _doSelect = doSelect;
        _nEstH = nEstH;
        _limitIRR = limitIRR;
        _gainMargin = gainMargin;
    }


    size_t inputBlockLength() @property
    {
        immutable a = _distorter.inputBlockLength,
                  b = (_nFFT + _nCP) * _nOS * (_doSelect ? _nEstH : 1);

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
    do{
        immutable size_t blk = this.inputBlockLength;
        immutable size_t dim = _distorter.outputDim;

        if(doLearning)
        {
            if(_doSelect && _importance !is null && _selectedBasisFuncs is null){
                _selectedBasisFuncs = new bool[][](_nFFT * _nOS, Dist.outputDim);
                _selectingIsSuccess = new bool[](_nFFT * _nOS);
                _selectingIsSuccess[] = true;
                _estimatedPower = new real[][](_nFFT * _nOS, Dist.outputDim);
                auto fftw = globalBankOf!(makeFFTWObject)[_nFFT * _nOS];

                C[][] freqX = new C[][](_nEstH, _nFFT * _nOS);
                C[][] freqY = new C[][](_nEstH, _nFFT * _nOS);
                foreach(idxOfEstH; 0 .. _nEstH){
                    fftw.inputs!R[] = input[idxOfEstH * (_nFFT + _nCP) * _nOS .. $][_nCP * _nOS .. (_nFFT + _nCP) * _nOS];
                    fftw.fft!R();
                    freqX[idxOfEstH][] = fftw.outputs!R[];

                    fftw.inputs!R[] = desired[idxOfEstH * (_nFFT + _nCP) * _nOS .. $][_nCP * _nOS .. (_nFFT + _nCP) * _nOS];
                    fftw.fft!R();
                    freqY[idxOfEstH][] = fftw.outputs!R[];
                }

                // MMSE推定する
                C[] freqH = new C[](_nFFT * _nOS);
                foreach(f; 0 .. _nFFT * _nOS){
                    real den = 0;
                    C num = C(0, 0);

                    foreach(idxOfEstH; 0 .. _nEstH){
                        den += freqX[idxOfEstH][f].sqAbs;
                        num += freqX[idxOfEstH][f].conj * freqY[idxOfEstH][f];
                    }

                    freqH[f] = cast(C)(num / den);
                }

                foreach(f; 0 .. _nFFT * _nOS){
                    R h = freqH[f].sqAbs;
                    R hinv = freqH[$-1-f].sqAbs;

                    foreach(p; 0 .. Dist.outputDim){
                        // 重要度と受信信号電力の積が，推定された干渉電力になる
                        immutable R ip = (){
                                real y = 0;
                                foreach(idxOfEstH; 0 .. _nEstH)
                                    y += freqY[idxOfEstH][f].sqAbs;
                                
                                return (y / _nEstH) * _importance[f][p];
                            }();
                        _estimatedPower[f][p] = ip;
                    }

                    // writefln("%s, %s, %s", f, h, hinv);

                    // この条件が満たされたとき，全基底関数を使用する
                    if(hinv > h * _limitIRR.gain^^2) {
                        _selectingIsSuccess[f] = false;
                        foreach(p; 0 .. Dist.outputDim)
                            _selectedBasisFuncs[f][p] = true;
                        
                        continue;
                    }

                    size_t cnt;
                    foreach(p; 0 .. Dist.outputDim) {
                        auto ip = _estimatedPower[f][p];

                        if(Dist.indexOfConjugated(p) >= p && !_selectedBasisFuncs[f][p]) {
                            // 干渉電力にマージンを設けて，ノイズ電力と比較する
                            if(ip * _gainMargin.gain^^2 > _noiseFloor){
                                _selectedBasisFuncs[f][p] = true;
                                _selectedBasisFuncs[f][Dist.indexOfConjugated(p)] = true;
                                cnt += 2;
                            }
                        }
                    }

                    if(cnt == 0){
                        _selectedBasisFuncs[f][0] = true;
                        cnt += 1;
                    }

                    // _states, _adaptersを更新する
                    _states[f] = MultiFIRState!C(cnt, 1);
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
                            _fftw.inputs!R[k] = bsym[k][p];

                        _fftw.fft!R();
                        _fftedBuffer[p][] = _fftw.outputs!R[];
                    }

                    // 受信信号の周波数成分
                    _fftw.inputs!R[] = dsym[];
                    _fftw.fft!R();
                    auto rxFFTed = _fftw.outputs!R();

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
                            _regenerator.frequencyResponse[f, p] = Complex!R(0, 0);
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
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
            _distorter.learn(signalGenerator(model).txBaseband);

        if(_doSelect)
            profiling(model, signalGenerator);
    }


    void profiling(M, Signals)(M originalModel, Signals delegate(M) genSignal)
    {
        import dffdd.filter.ls : LSAdapter, makeLSAdapter;

        // testcase == 0 : 重要度を計算する
        // testcase == 1 : デバッグ用に_actualPowerを計算する
        foreach(testcase; [1, 0])
        {
            auto model = originalModel;

            if(testcase == 0)
            {
                // 同軸線路を使用するように設定する
                model.useCoaxialCableAsChannel();
                model.rndSeed += 114514;
            }

            // 最初にノイズ電力の計測
            // _noiseFloor = signals.noise.save.map!`a.re^^2 + a.im^^2`.take(1024*1024).sum() / (1024 * 1024);
            {
                auto signals = genSignal(model);

                auto n = signals.noise.save;
                auto fftw = globalBankOf!(makeFFTWObject)[_nFFT * _nOS];
                R[] buf = new R[_nFFT * _nOS];
            
                foreach(ref e; buf)
                    e = 0;

                // auto buf = new Complex!float[_nOS * _nFFT];
                foreach(i; 0 .. 1024){
                    foreach(j; 0 .. _nOS * _nFFT){
                        fftw.inputs!R[j] = cast(Complex!R)n.front;
                        n.popFront();
                    }

                    fftw.fft!R();
                    foreach(j; 0 .. _nOS * _nFFT)
                        buf[j] += fftw.outputs!R[j].sqAbs;
                }

                foreach(j; 0 .. _nOS * _nFFT)
                    buf[j] /= 1024;

                _noiseFloor = buf.sum() / (_nFFT * _nOS);
            }

            if(testcase == 0)
            {
                // LNAダイナミックレンジ最大まで使用する
                model.INR = model.lna.DR;
            }

            // {
            //     auto model4IRR = model;
            //     model.ofdm.numOfSubcarrier = 1;
            //     model.noImpairementsOnTX();
            // }

            auto signals = genSignal(model);
            signals.useSWPOFDM = true;
            signals.ignoreDesired = true;

            // 次のようなフィルタを学習してみる
            // + 適応アルゴリズム : 最小二乗法
            // + 使用シンボル数：1024シンボル
            auto testfilter = new SimpleFrequencyDomainParallelHammersteinFilter!(C, Dist, LSAdapter!(MultiFIRState!C))
                            (   _distorter,
                                (size_t i, bool b, MultiFIRState!C s) => makeLSAdapter(s, 1024),
                                _scMap, _nFFT, _nCP, _nOS, _sampFreq, false);

            immutable spb = this.inputBlockLength / ((_nFFT + _nCP) * _nOS);    // 1ブロックあたりのシンボル数

            C[] testTX = new C[this.inputBlockLength],
                testRX = new C[this.inputBlockLength],
                testER = new C[this.inputBlockLength];

            // 1024シンボル使用して学習する
            foreach(i; 0 .. 1024 / spb){
                signals.fillBuffer!(["txBaseband", "received"])(testTX, testRX);
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
            enum size_t Navg = 1024;
            foreach(i; 0 .. Navg){
                signals.fillBuffer!(["txBaseband", "received"])(testTX, testRX);
                _distorter(testTX, disted);

                // 非線形送信信号を周波数領域に変換
                foreach(p; 0 .. _distorter.outputDim){
                    foreach(k; 0 .. _nFFT * _nOS)
                        _fftw.inputs!R[k] = disted[_nCP * _nOS + k][p];

                    _fftw.fft!R();

                    foreach(k; 0 .. _nFFT * _nOS){
                        if(i == 0) firstPowerOfSI[k][p] = _fftw.outputs!R[k].sqAbs;
                        disted[k][p] = _fftw.outputs!R[k] * weight[k][p];
                    }
                }

                // 受信信号の周波数成分
                _fftw.inputs!R[] = testRX[_nCP * _nOS .. $];
                _fftw.fft!R();
                auto rxFFTed = _fftw.outputs!R();

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

            if(testcase == 0){          // 重要度を計算する場合なら
                _importance = powerOfSI;
                foreach(f; 0 .. _nFFT * _nOS) foreach(p; 0 .. _distorter.outputDim){
                    _importance[f][p] /= powerOfRCV[f];
                }
            }else if(testcase == 1){    // _actualPowerを算出するだけなら
                _actualPower = powerOfSI;
                foreach(ref es; _actualPower) foreach(ref e; es) e /= Navg;
                _requiredBasisFuncs = new bool[][](_nFFT * _nOS, _distorter.outputDim);
                foreach(f; 0 .. _nFFT * _nOS) foreach(p; 0 .. _distorter.outputDim)
                    _requiredBasisFuncs[f][p] = (_actualPower[f][p] >= _noiseFloor);
            }
        }
    }


    void saveInfoToDir(string dir)
    {
        import std.file : mkdirRecurse;
        import std.stdio;
        import std.path;
        import std.json;
        import std.conv;

        dir = buildPath(dir, "filterSpec");
        mkdirRecurse(dir);

        if(_doSelect){
            auto file = File(buildPath(dir, "significance.csv"), "w");
            // 重要度
            foreach(f; 0 .. _nFFT * _nOS)
                file.writefln("%s, %(%s, %)", f, _importance[f].map!(a => 10*log10(a)));

            // サブキャリアの使用状況
            file = File(buildPath(dir, "basisfuncs_per_subcarrier.csv"), "w");
            foreach(f; 0 .. _nFFT * _nOS)
                file.writefln("%s, %(%s, %)", f, _selectedBasisFuncs[f].map!`a ? 1 : 0`.zip(iota(1, _selectedBasisFuncs[f].length+1)).map!`a[0]*a[1]`);

            // selectingIsSuccess
            file = File(buildPath(dir, "success_estimation_for_selecting.csv"), "w");
            foreach(f; 0 .. _nFFT * _nOS)
                file.writefln("%s, %s", f, _selectingIsSuccess[f] ? 1 : 0);

            // estimatedPower
            file = File(buildPath(dir, "estimated_power.csv"), "w");
            foreach(f; 0 .. _nFFT * _nOS)
                file.writefln("%s, %(%s, %)", f, _estimatedPower[f].map!(a => a <= 0 ? -400 : 10*log10(a)));

            // actualPower
            file = File(buildPath(dir, "actual_power.csv"), "w");
            foreach(f; 0 .. _nFFT * _nOS)
                file.writefln("%s, %(%s, %)", f, _actualPower[f].map!(a => a <= 0 ? -400 : 10*log10(a)));

            file = File(buildPath(dir, "actual_required_basis_funcs.csv"), "w");
            foreach(f; 0 .. _nFFT * _nOS)
                file.writefln("%s, %(%s, %)", f, _requiredBasisFuncs[f].map!`a ? 1 : 0`.zip(iota(1, _requiredBasisFuncs[f].length+1)).map!`a[0]*a[1]`);
        }

        auto jv = this.info;
        auto file = File(buildPath(dir, "info.json"), "w");
        file.writeln(jv.toPrettyString());
    }


    JSONValue info()
    {
        import std.conv;

        JSONValue jv = [ "type": typeof(this).stringof ];
        if(_doSelect){
            size_t[] bfcounter(bool[][] arr) {
                size_t[] cnts = new size_t[_distorter.outputDim+1];
                foreach(f; 0 .. _nFFT * _nOS){
                    auto cnt = arr[f].map!`a ? 1 : 0`.sum();
                    cnts[cnt] += 1;
                }
                return cnts;
            }

            jv["limitIRR"] = _limitIRR.to!string;
            jv["gainMargin"] = _gainMargin.to!string;
            jv["numOfEstimationOfH"] = _nEstH;
            jv["selectedBasisFuncs"] = _selectedBasisFuncs;
            jv["selectedBasisFuncsCounts"] = bfcounter(_selectedBasisFuncs);
            jv["noiseFloor"] = (10*log10(_noiseFloor)).dB.to!string;
            jv["selectingIsSuccess"] = _selectingIsSuccess;
            jv["estimatedPower"] = _estimatedPower;
            jv["actualPower"] = _actualPower;
            jv["actualRequiredBasisFuncs"] = _requiredBasisFuncs;
            jv["actualRequiredBasisFuncsCounts"] = bfcounter(_requiredBasisFuncs);
        }
        jv["scMap"] = _scMap.map!(a => a ? 1 : 0).array();
        jv["nFFT"] = _nFFT;
        jv["nCP"] = _nCP;
        jv["nOS"] = _nOS;
        jv["sampFreq"] = _sampFreq;
        jv["doSelect"] = _doSelect;

        return jv;
    }


  private:
    FFTWObject!Complex _fftw;
    Dist _distorter;
    OverlapSaveRegenerator!C _regenerator;
    Adapter delegate(size_t, bool, MultiFIRState!C) _genAdapter;
    immutable(bool[]) _scMap;
    immutable size_t _nFFT, _nCP, _nOS;
    immutable real _sampFreq;
    immutable bool _doSelect;
    immutable size_t _nEstH;
    immutable Gain _limitIRR, _gainMargin;
    MultiFIRState!C[] _states;
    Adapter[] _adapters;
    C[][] _distortedBuffer;
    C[][] _fftedBuffer;

    // 重要度計算
    real[][] _importance;
    bool[][] _selectedBasisFuncs;
    R _noiseFloor = 0;
    // 重要度でのログ
    bool[] _selectingIsSuccess;     // 各周波数で，重要度が指標として使用可能ならtrue, それ以外ならfalse
    real[][] _estimatedPower;       // 重要度から推定された各基底の電力
    real[][] _actualPower;
    bool[][] _requiredBasisFuncs;
}


final class SimpleFrequencyDomainDCMHammersteinFilterType2(C, Dist, alias genAdaptor, Flag!"isParallel" isParallel = No.isParallel)
if(isBlockConverter!(Dist, C, C[]))
{
    alias R = typeof(C.init.re);
    enum bool doesDistHaveLearn = is(typeof((Dist dist, C[] r){ dist.learn(r); }));

    import dffdd.filter.regenerator : OverlapSaveRegenerator;

    this(Dist dist, in bool[] subcarrierMap, size_t nFFT, size_t nCP, size_t nOS)
    {
        immutable dim = dist.outputDim;

        _fftw = makeFFTWObject!Complex(nFFT * nOS);
        _distorter = dist;
        auto sliceState = new C[nFFT * nOS * dim].sliced(nFFT * nOS, dim);
        sliceState[] = C(0);

        _regenerator = new OverlapSaveRegenerator!C(sliceState);

        foreach(i; 0 .. nFFT * nOS){
            MultiFIRState!C[] sts;
            typeof(_adapters[0]) ads;
            foreach(p; 0 .. dim){
                auto state = MultiFIRState!C(1, 1);
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
    do{
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
                            _fftw.inputs!R[k] = bsym[k][p];

                        _fftw.fft!R();
                        _fftedBuffer[p][] = _fftw.outputs!R[];
                    }

                    // 受信信号の周波数成分
                    _fftw.inputs!R[] = dsym[];
                    _fftw.fft!R();
                    auto rxFFTed = _fftw.outputs!R();

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
    MultiFIRState!C[][] _states;
    typeof(genAdaptor(0, true, 0, _states[0][0]))[][] _adapters;
    C[][] _distortedBuffer;
    C[][] _fftedBuffer;
}


final class SimpleFrequencyDomainDCMHammersteinFilterType1(C, Dist, alias genAdaptor, Flag!"isParallel" isParallel = No.isParallel)
{
    alias R = typeof(C.init.re);
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
    do{
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
