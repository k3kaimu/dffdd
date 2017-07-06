module dffdd.filter.freqdomain;

import dffdd.filter.primitives;
import dffdd.filter.ls;
import dffdd.filter.state;
import dffdd.filter.traits;
import dffdd.filter.filter;
import dffdd.utils.unit;
import dffdd.utils.fft;

import std.algorithm;
import std.complex;
import std.json;
import std.math;
import std.numeric;
import std.range;
import std.typecons;

import carbon.math : complexZero;


enum bool isFrequencyDomainMISOStateAdapter(T, C = Complex!float) = is(typeof((T adapter){
    C[][] distX;    // 送信信号とその歪信号の各周波数成分，各周波数で必要な次元数にすでに減らされている
    C[] desired;    // 所望信号の各周波数成分
    size_t[] dims;  // 各周波数成分は何次元か

    adapter.setDimensionForEachFreq(dims);
    adapter.adapt(distX, desired);
    
    C[] output;     // 各周波数での推定値
    adapter.regenerate(distX, output);
}));


final class FrequencyDomainParallelHammersteinStateAdapter(C, Adapter)
{
    this(Adapter delegate(size_t i, bool b, MultiFIRState!C) genAdapter, in bool[] subcarrierMap, size_t maxDim)
    {
        _nFFT = subcarrierMap.length;
        _scMap = subcarrierMap.dup;
        _genAdapter = genAdapter;

        foreach(i; 0 .. _nFFT){
            auto state = MultiFIRState!(Complex!float)(maxDim, 1);
            _states ~= state;
            _adapters ~= _genAdapter(i, _scMap[i], state);
        }
    }


    void setDimensionForEachFreq(size_t[] dims)
    {
        _states.length = 0;
        _adapters.length = 0;

        foreach(i; 0 .. _nFFT){
            immutable dim = dims[i];

            auto state = MultiFIRState!(Complex!float)(dim, 1);
            _states ~= state;
            _adapters ~= _genAdapter(i, _scMap[i], state);
        }
    }


    void adapt(C[][] distX, C[] desired)
    {
        foreach(f; 0 .. _nFFT){
            assert(distX[f].length == _states[f].numOfFIR);
            if(distX[f].length == 0) continue;

            _states[f].update(distX[f]);
            C error = desired[f] - _states[f].output;
            _adapters[f].adapt(_states[f], error);
        }
    }


    void regenerate(C[][] distX, ref C[] dst)
    {
        if(dst.length != _nFFT) dst.length = _nFFT;

        foreach(f; 0 .. _nFFT){
            dst[f] = complexZero!C;

            foreach(i; 0 .. _states[f].numOfFIR)
                dst[f] += _states[f].weight[0][i] * distX[f][i];
        }
    }


  private:
    size_t _nFFT;   // nOS * nFFTの値
    immutable(bool[]) _scMap;
    Adapter delegate(size_t i, bool b, MultiFIRState!C) _genAdapter;
    MultiFIRState!C[] _states;
    Adapter[] _adapters;
}


final class FrequencyDomainDCMHammersteinStateAdapter(C, Adapter, Flag!"isParallel" isParallel = No.isParallel)
{
    this(Adapter delegate(size_t i, bool b, size_t p, MultiFIRState!C) genAdapter, in bool[] subcarrierMap, size_t maxDim)
    {
        _nFFT = subcarrierMap.length;
        _scMap = subcarrierMap.dup;
        _genAdapter = genAdapter;

        foreach(i; 0 .. _nFFT){
            MultiFIRState!C[] sts;
            Adapter[] ads;
            foreach(p; 0 .. maxDim){
                auto state = MultiFIRState!(Complex!float)(1, 1);
                sts ~= state;
                ads ~= _genAdapter(i, _scMap[i], p, state);
            }

            _states ~= sts;
            _adapters ~= ads;
        }
    }


    void setDimensionForEachFreq(size_t[] dims)
    {
        _states.length = 0;
        _adapters.length = 0;

        foreach(i; 0 .. _nFFT){
            immutable dim = dims[i];

            MultiFIRState!C[] sts;
            Adapter[] ads;
            foreach(p; 0 .. dim){
                auto state = MultiFIRState!(Complex!float)(1, 1);
                sts ~= state;
                ads ~= _genAdapter(i, _scMap[i], p, state);
            }

            _states ~= sts;
            _adapters ~= ads;
        }
    }


    void adapt(C[][] distX, C[] desired)
    {
        foreach(f; 0 .. _nFFT){
            immutable dim = distX[f].length;

            C error = desired[f];
            foreach(p; 0 .. dim){
                _states[f][p].update(distX[f][p]);
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


    void regenerate(C[][] distX, ref C[] dst)
    {
        if(dst.length != _nFFT) dst.length = _nFFT;

        foreach(f; 0 .. _nFFT){
            C e = complexZero!C;

            foreach(p, ref st; _states[f])
                e += st.weight[0, 0] * distX[f][p];

            dst[f] = e;
        }
    }


  private:
    size_t _nFFT;   // nOS * nFFTの値
    immutable(bool[]) _scMap;
    Adapter delegate(size_t i, bool b, size_t p, MultiFIRState!C) _genAdapter;
    MultiFIRState!C[][] _states;
    Adapter[][] _adapters;
}


final class OverlapSaveRegenerator2(C)
{
    this(size_t maxDim, size_t nFFT)
    {
        _fftw = makeFFTWObject!Complex(nFFT);
        _nFFT = nFFT;
        _maxDim = maxDim;
        _buffer = new C[](nFFT);
        _inputs = new C[][](maxDim, nFFT);
        foreach(ref es; _inputs) foreach(ref e; es) e = complexZero!C;
        _distFFTBuf = new C[][](maxDim, nFFT);
        _distFFTBufSelected = new C[][](nFFT, maxDim);
    }


    void setSelectedForEachFreq(in bool[][] selected)
    {
        _selected.length = 0;
        foreach(f, e; selected){
            _selected ~= e.dup;

            size_t cnt = e.map!(a => a ? 1 : 0).sum();
            _distFFTBufSelected[f].length = cnt;
        }
    }


    void apply(MISOState)(ref MISOState state, const(C[])[] tx, ref C[] output)
    in{
        foreach(es; tx) foreach(e; es){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    out{
        foreach(e; output){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    body{
        immutable batchLen = _nFFT / 2;

        output.length = tx.length;
        size_t remain = tx.length;

        auto ops = output;
        while(remain){
            auto size = min(remain, batchLen);
            applyImpl(state, tx[0 .. size], ops[0 .. size]);
            tx = tx[size .. $];
            ops = ops[size .. $];
            remain -= size;
        }
    }


    void applyImpl(MISOState)(ref MISOState state, in C[][] tx, C[] output)
    in{
        assert(tx.length <= _nFFT);
        assert(tx.length == output.length);
        foreach(es; _inputs) foreach(e; es){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }

        foreach(es; tx) foreach(e; es){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    out{
        foreach(e; output){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    body{
        immutable size_t size = tx.length;
        output[] = complexZero!C;

        // _bufferをゼロ初期化
        foreach(ref e; _buffer) e = complexZero!C;

        foreach(p; 0 .. _maxDim)
        {
            // _inputsをシフトして，後ろにsize分の空きを作る
            foreach(i; 0 .. _inputs[p].length - size)
                _inputs[p][i] = _inputs[p][i + size];

            // 後ろにできた空き部分にコピーする
            foreach(i, ref e; _inputs[p][$ - size .. $])
                e = tx[i][p];

            // 信号xを周波数領域に変換しXを得る
            auto ips = _fftw.inputs!float;
            auto ops = _fftw.outputs!float;
            ips[] = _inputs[p][];
            _fftw.fft!float();

            _distFFTBuf[p][] = ops[];
        }

        // _selectedで指定された基底だけ抽出する
        foreach(f; 0 .. _nFFT){
            size_t cnt;
            foreach(p; 0 .. _maxDim){
                if(_selected is null || _selected[f][p]){
                    _distFFTBufSelected[f][cnt] = _distFFTBuf[p][f];
                    ++cnt;
                }
            }
        }

        // XとHの積を計算する
        auto ips = _fftw.inputs!float;
        state.regenerate(_distFFTBufSelected, ips);
        assert(ips.ptr == _fftw.inputs!float.ptr);

        // y = IFFT{XH}の計算
        // IFFT
        // _fftw.inputs!float[] = _buffer[];
        _fftw.ifft!float();
        // 結果の後ろをoutputに足す
        output[] += _fftw.outputs!float[$ - size .. $];
    }


  private:
    FFTWObject!Complex _fftw;
    size_t _nFFT, _maxDim;
    C[][] _inputs;
    C[] _buffer;
    immutable(bool[])[] _selected;
    C[][] _distFFTBuf, _distFFTBufSelected;
}


final class FrequencyDomainHammersteinFilter(C, Dist, StateAdapter)
if(isBlockConverter!(Dist, C, C[]) && isFrequencyDomainMISOStateAdapter!(StateAdapter, C))
{
    enum bool doesDistHaveLearn = is(typeof((Dist dist, C[] r){ dist.learn(r); }));

    this(Dist dist, StateAdapter stateAdapter, in bool[] subcarrierMap, size_t nFFT, size_t nCP, size_t nOS, real sampFreq,
        bool doSelect = false, bool doComplexSelect = false, size_t nEstH = 2, Gain limitIRR = 20.dB, Gain gainMargin = 6.dB)
    {
        immutable dim = dist.outputDim;

        _fftw = makeFFTWObject!Complex(nFFT * nOS);
        _distorter = dist;
        _stateAdapter = stateAdapter;
        _regenerator = new OverlapSaveRegenerator2!C(dim, nFFT * nOS);

        _scMap = subcarrierMap.dup;
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;
        _nSC = subcarrierMap.map!(a => a ? 1 : 0).sum();
        _sampFreq = sampFreq;
        _distortedBuffer = new C[][](this.inputBlockLength, dim);
        _fftedBuffer = new C[][](dim, _nFFT * _nOS);
        _doSelect = doSelect;
        _doComplexSelect = doComplexSelect;
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
    body{
        immutable size_t blk = this.inputBlockLength;
        immutable size_t dim = _distorter.outputDim;

        if(doLearning)
        {
            if(_doSelect && _importance !is null && _selectedBasisFuncs is null){
                _selectedBasisFuncs = new bool[][](_nFFT * _nOS, Dist.outputDim);
                _selectedBasisFuncsDim = new size_t[_nFFT * _nOS];
                _selectingIsSuccess = new bool[](_nFFT * _nOS);
                _selectingIsSuccess[] = true;
                _estimatedPower = new real[][](_nFFT * _nOS, Dist.outputDim);
                //auto fftw = globalBankOf!(makeFFTWObject)[_nFFT * _nOS];

                C[][] freqX = new C[][](_nEstH, _nFFT * _nOS);
                C[][] freqY = new C[][](_nEstH, _nFFT * _nOS);
                foreach(idxOfEstH; 0 .. _nEstH){
                    _fftw.inputs!float[] = input[idxOfEstH * (_nFFT + _nCP) * _nOS .. $][_nCP * _nOS .. (_nFFT + _nCP) * _nOS];
                    _fftw.fft!float();
                    freqX[idxOfEstH][] = _fftw.outputs!float[];

                    _fftw.inputs!float[] = desired[idxOfEstH * (_nFFT + _nCP) * _nOS .. $][_nCP * _nOS .. (_nFFT + _nCP) * _nOS];
                    _fftw.fft!float();
                    freqY[idxOfEstH][] = _fftw.outputs!float[];
                }

                if(_doComplexSelect)
                {
                    foreach(f; 0 .. _nFFT * _nOS){
                        size_t cnt;
                        foreach(p; 0 .. Dist.outputDim){
                            // 受信信号電力E[|Y[f]|^2]
                            immutable float pY0 = (){
                                real y = 0;
                                    foreach(idxOfEstH; 0 .. _nEstH)
                                        y += freqY[idxOfEstH][f].sqAbs;
                                
                                return y / _nEstH;
                            }();

                            // 受信信号電力E[|Y[-f]|^2]
                            immutable float pY1 = (){
                                real y = 0;
                                    foreach(idxOfEstH; 0 .. _nEstH)
                                        y += freqY[idxOfEstH][f == 0 ? 0 : $-f].sqAbs;
                                
                                return y / _nEstH;
                            }();

                            immutable iqcoef = (1 / _limitIRR.gain)^^2;

                            immutable float pY2 = (pY0 + iqcoef * pY1 + 2 * sqrt(pY0 * iqcoef * pY1)) / (1 - iqcoef)^^2;
                            immutable float pY3 = (pY1 + iqcoef * pY0 + 2 * sqrt(pY1 * iqcoef * pY0)) / (1 - iqcoef)^^2;

                            // 重要度と受信信号電力の積が，推定された干渉電力になる
                            immutable float ipR = (){
                                    return pY2 * _importance[f][p];
                            }();

                            immutable float ipI = (){
                                    return pY3 * _importance[f == 0 ? 0 : $-f][_distorter.indexOfConjugated(p)]
                                               * iqcoef;
                            }();

                            _estimatedPower[f][p] = ipR + ipI + 2 * sqrt(ipR * ipI);

                            if(_estimatedPower[f][p] * _gainMargin.gain^^2 > _noiseFloor){
                                _selectedBasisFuncs[f][p] = true;
                                cnt += 1;
                            }
                        }

                        _selectedBasisFuncsDim[f] = cnt;
                    }
                }
                else
                {
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
                        float h = freqH[f].sqAbs;
                        float hinv = freqH[$-1-f].sqAbs;

                        foreach(p; 0 .. Dist.outputDim){
                            // 重要度と受信信号電力の積が，推定された干渉電力になる
                            immutable float ip = (){
                                    real y = 0;
                                    foreach(idxOfEstH; 0 .. _nEstH)
                                        y += freqY[idxOfEstH][f].sqAbs;
                                    
                                    return (y / _nEstH) * _importance[f][p];
                                }();
                            _estimatedPower[f][p] = ip;
                        }

                        // この条件が満たされたとき，全基底関数を使用する
                        size_t cnt;
                        if(hinv > h * _limitIRR.gain^^2) {
                            cnt = Dist.outputDim;
                            _selectingIsSuccess[f] = false;
                            foreach(p; 0 .. Dist.outputDim)
                                _selectedBasisFuncs[f][p] = true;
                            
                            goto LnextFreq;
                        }

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

                      LnextFreq:
                        if(cnt == 0){
                            _selectedBasisFuncs[f][0] = true;
                            cnt += 1;
                        }

                        assert(cnt != 0);

                        // _states, _adaptersを更新する
                        // _states[f] = MultiFIRState!(Complex!float)(cnt, 1);
                        // _adapters[f] = _genAdapter(f, _scMap[f], _states[f]);
                        _selectedBasisFuncsDim[f] = cnt;
                    }
                }

                _stateAdapter.setDimensionForEachFreq(_selectedBasisFuncsDim);
                _regenerator.setSelectedForEachFreq(_selectedBasisFuncs);
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

                    C[][] distX = new C[][](_nFFT * _nOS);
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

                        distX[f] = vX2.dup;
                    }

                    _stateAdapter.adapt(distX, rxFFTed);
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

            _regenerator.apply(_stateAdapter, _distortedBuffer, ers);
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
            import std.stdio;
            auto model = originalModel;

            with(model){
                useDesiredBaseband = false;
                useTxBaseband = false;
                useTxBasebandSWP = true;        // 使う
                useDesiredPADirect = false;
                usePADirect = false;
                usePADirectSWP = false;
                useReceivedDesired = false;
                useReceivedSI = false;
                useReceivedSISWP = true;        // 使う
                useReceived = false;
                useNoise = true;                // 使う
            }

            if(testcase == 0)
            {
                // 同軸線路を使用するように設定する
                model.useCoaxialCableAsChannel();
                model.rndSeed += 114514;
            }

            // 最初にノイズ電力の計測
            {
                model.useTxBasebandSWP = false;
                model.useReceivedSISWP = false;
                auto signals = genSignal(model);
                model.useTxBasebandSWP = true;
                model.useReceivedSISWP = true;

                auto n = signals.noise.save;
                // auto fftw = globalBankOf!(makeFFTWObject)[_nFFT * _nOS];
                float[] buf = new float[_nFFT * _nOS];
            
                foreach(ref e; buf)
                    e = 0;

                // auto buf = new Complex!float[_nOS * _nFFT];
                immutable size_t nSize = 64;
                foreach(i; 0 .. nSize){
                    foreach(j; 0 .. _nOS * _nFFT){
                        _fftw.inputs!float[j] = cast(Complex!float)n.front;
                        n.popFront();
                    }

                    _fftw.fft!float();
                    foreach(j; 0 .. _nOS * _nFFT)
                        buf[j] += _fftw.outputs!float[j].sqAbs;
                }

                foreach(j; 0 .. _nOS * _nFFT)
                    buf[j] /= nSize;

                _noiseFloor = buf.sum() / (_nFFT * _nOS);
            }

            if(testcase == 0)
            {
                // LNAダイナミックレンジ最大まで使用する
                model.INR = model.lna.DR;
            }

            C[][] weight;
            Tuple!(real[][], real[]) expectedPowers;
            real[][] psiPower;
            real[] gnl = new real[_distorter.outputDim];
            auto signals = genSignal(model);
            if(_doComplexSelect && testcase == 0)
            {
                _iqRX = estimateIQCoefs(estimateCFR(signals.save(), C(0, 0)))[1];
                weight = estimateCFR(signals.save(), _iqRX);
                foreach(p; 0 .. _distorter.outputDim)
                {
                    gnl[p] = 0;
                    foreach(f; 0 .. _nFFT * _nOS)
                        if(_scMap[f]) gnl[p] += weight[f][p].sqAbs / weight[f][0].sqAbs;
                    
                    gnl[p] = sqrt(gnl[p] / _nSC);
                }

                //expectedPowers = calculateExpectedPower(signals, weight, C(0, 0));
                psiPower = calculateExpectedPower(signals.save(), weight, C(0, 0))[2];
            }
            else
            {
                // H_{p,q}[f] を取得する
                weight = estimateCFR(signals.save(), C(0, 0));
                expectedPowers[0 .. 2] = calculateExpectedPower(signals.save(), weight, C(0, 0))[0 .. 2];
            }

            if(testcase == 0){          // 重要度を計算する場合なら
                if(_doComplexSelect){
                    // _importance = psiPower.dup;
                    _importance = new real[][](_nFFT * _nOS, _distorter.outputDim);
                    foreach(f; 0 .. _nFFT * _nOS) foreach(p; 0 .. _distorter.outputDim){
                        _importance[f][p] = psiPower[f][p] / psiPower[f][0];
                        _importance[f][p] *= gnl[p]^^2;
                    }
                }else{
                    _importance = new real[][](_nFFT * _nOS, _distorter.outputDim);
                    foreach(f; 0 .. _nFFT * _nOS) foreach(p; 0 .. _distorter.outputDim){
                        _importance[f][p] = expectedPowers[0][f][p] / expectedPowers[0][f][0];
                    }
                }
            }else if(testcase == 1){    // _actualPowerを算出するだけなら
                _actualPower = expectedPowers[0];
                // foreach(ref es; _actualPower) foreach(ref e; es) e /= Navg;
                _requiredBasisFuncs = new bool[][](_nFFT * _nOS, _distorter.outputDim);
                foreach(f; 0 .. _nFFT * _nOS) foreach(p; 0 .. _distorter.outputDim)
                    _requiredBasisFuncs[f][p] = (_actualPower[f][p] >= _noiseFloor);
            }
        }
    }


    C[][] estimateCFR(Signals)(Signals signals, C iqRX = C(0, 0))
    {
        // 次のようなフィルタを学習してみる
        // + 適応アルゴリズム : 最小二乗法
        // + 使用シンボル数：1024シンボル
        auto testfilter = new FrequencyDomainHammersteinFilter!(Complex!float, Dist, FrequencyDomainParallelHammersteinStateAdapter!(Complex!float, LSAdapter!(MultiFIRState!(Complex!float))))
                            (_distorter,
                            new FrequencyDomainParallelHammersteinStateAdapter!(Complex!float, LSAdapter!(MultiFIRState!(Complex!float)))
                                ((size_t i, bool b, MultiFIRState!C s) => makeLSAdapter(s, 64), _scMap, _distorter.outputDim)
                            ,
                            _scMap, _nFFT, _nCP, _nOS, _sampFreq, false, false);

        immutable spb = this.inputBlockLength / ((_nFFT + _nCP) * _nOS);    // 1ブロックあたりのシンボル数

        C[] testTX = new C[this.inputBlockLength],
            testRX = new C[this.inputBlockLength],
            testER = new C[this.inputBlockLength];

        // 64シンボル使用して学習する
        foreach(i; 0 .. 64 / spb){
            signals.fillBuffer!(["txBasebandSWP", "receivedSISWP"])(testTX, testRX);

            testfilter.apply!(Yes.learning)(testTX, testRX, testER);
        }

        // 各周波数で，学習完了後の重みを得る
        C[][] weight = new C[][](_nFFT * _nOS, _distorter.outputDim);
        foreach(f; 0 .. _nFFT * _nOS)
            foreach(p; 0 .. _distorter.outputDim){
                immutable w1 = testfilter._stateAdapter._states[f].weight[0][p];
                immutable w2 = testfilter._stateAdapter._states[$-1-f].weight[0][_distorter.indexOfConjugated(p)];

                weight[f][p] = (w1 - iqRX * w2.conj) / (1 - iqRX.sqAbs);
            }

        return weight;
    }


    C[2] estimateIQCoefs(C[][] weights)
    {
        import std.experimental.ndslice;
        import dffdd.utils.linalg;

        Slice!(2, C*) mx = new C[2 * _nSC].sliced(2, _nSC);
        C[] y = new C[_nSC];
        size_t cnt;
        foreach(f; 0 .. _nFFT * _nOS) if(_scMap[f]) {
            mx[0, cnt] = weights[f][0];
            mx[1, cnt] = weights[f == 0 ? 0 : $-f][0].conj;
            y[cnt] = weights[f][1];
            ++cnt;
        }

        auto h = leastSquareEstimate(mx, y);
        return [h[0], h[1]];
    }


    Tuple!(real[][], real[], real[][]) calculateExpectedPower(Signals)(Signals signals, in C[][] weight, C iqRX = C(0, 0))
    {
        C[] testTX = new C[this.inputBlockLength],
            testRX = new C[this.inputBlockLength],
            testER = new C[this.inputBlockLength];

        real[][] powerOfPsi = new real[][](_nFFT * _nOS, _distorter.outputDim);
        real[][] powerOfSI = new real[][](_nFFT * _nOS, _distorter.outputDim);
        real[] powerOfRCV = new real[](_nFFT * _nOS);

        foreach(f; 0 .. _nFFT * _nOS){
            powerOfRCV[f] = 0;
            foreach(p; 0 .. _distorter.outputDim){
                powerOfSI[f][p] = 0;
                powerOfPsi[f][p] = 0;
            }
        }

        testTX = testTX[0 .. (_nFFT + _nCP) * _nOS];
        testRX = testRX[0 .. (_nFFT + _nCP) * _nOS];
        auto disted = new C[][]((_nFFT + _nCP) * _nOS, _distorter.outputDim);

        // 重要度を計算する
        enum size_t Navg = 64;
        foreach(i; 0 .. Navg){
            signals.fillBuffer!(["txBasebandSWP", "receivedSISWP"])(testTX, testRX);
            _distorter(testTX, disted);

            // 非線形送信信号を周波数領域に変換
            foreach(p; 0 .. _distorter.outputDim){
                foreach(k; 0 .. _nFFT * _nOS)
                    _fftw.inputs!float[k] = disted[_nCP * _nOS + k][p];

                _fftw.fft!float();

                foreach(k; 0 .. _nFFT * _nOS){
                    // if(i == 0) firstPowerOfSI[k][p] = _fftw.outputs!float[k].sqAbs;
                    powerOfPsi[k][p] += _fftw.outputs!float[k].sqAbs;
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

            // if(i == 0)
            //     firstPowerOfRCV = powerOfRCV.dup;
        }

        foreach(f; 0 .. _nFFT * _nOS){
            foreach(p; 0 .. _distorter.outputDim){
                powerOfSI[f][p] /= Navg;
                powerOfPsi[f][p] /= Navg;
            }

            powerOfRCV[f] /= Navg;
        }

        return tuple(powerOfSI, powerOfRCV, powerOfPsi);
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
    StateAdapter _stateAdapter;
    OverlapSaveRegenerator2!C _regenerator;
    immutable(bool[]) _scMap;
    immutable size_t _nFFT, _nCP, _nOS, _nSC;
    immutable real _sampFreq;
    immutable bool _doSelect, _doComplexSelect;
    immutable size_t _nEstH;
    immutable Gain _limitIRR, _gainMargin;
    C _iqRX;
    C[][] _distortedBuffer;
    C[][] _fftedBuffer;

    // 重要度計算
    real[][] _importance;
    bool[][] _selectedBasisFuncs;
    size_t[] _selectedBasisFuncsDim;
    float _noiseFloor = 0;
    // 重要度でのログ
    bool[] _selectingIsSuccess;     // 各周波数で，重要度が指標として使用可能ならtrue, それ以外ならfalse
    real[][] _estimatedPower;       // 重要度から推定された各基底の電力
    real[][] _actualPower;
    bool[][] _requiredBasisFuncs;
}

unittest
{
    alias test1 = FrequencyDomainHammersteinFilter!(Complex!float, PADistorter!(Complex!float, 3), FrequencyDomainParallelHammersteinStateAdapter!(Complex!float, LSAdapter!(MultiFIRState!(Complex!float))));
    alias test2 = FrequencyDomainHammersteinFilter!(Complex!float, PADistorter!(Complex!float, 3), FrequencyDomainDCMHammersteinStateAdapter!(Complex!float, LSAdapter!(MultiFIRState!(Complex!float))));

    if(0){
        test1 obj;
        obj.apply!(Yes.learning)(null, null,  null);
    }
}


final class IQInversionSuccessiveInterferenceCanceller(C, size_t P)
{
    import std.experimental.ndslice;
    import dffdd.utils.linalg;


    this(size_t nWLLearning, size_t nIter, in bool[] subcarrierMap, size_t nFFT, size_t nCP, size_t nOS)
    {
        _scMap = subcarrierMap.dup;
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;
        _nWLLearning = nWLLearning;
        _nIter = nIter;

        _fftw = makeFFTWObject!Complex(_nFFT * _nOS);
        _regenerator = new OverlapSaveRegenerator2!C(1, _nFFT * _nOS);

        _iqTX = 0;
        _iqRX = 0;
        _paCoefs[0] = 1;
        foreach(ref e; _paCoefs[1 .. $]) e = 0;
        _channelFreqResponse = new Complex!real[_nFFT * _nCP];
        foreach(ref e; _channelFreqResponse) e = 0;

        _buffer4Regen = new C[(_nFFT + _nCP) * _nOS];
    }


    size_t inputBlockLength() @property
    {
        return (_nFFT + _nCP) * _nOS;
    }


    void apply(Flag!"learning" doLearning)(in C[] input, in C[] desired, C[] errors)
    {
        immutable nSym = (_nFFT + _nCP) * _nOS;

        foreach(i; 0 .. input.length / nSym)
        {
            applyForEachSymbol!doLearning(input[i*nSym .. (i+1)*nSym], desired[i*nSym .. (i+1)*nSym], errors[i*nSym .. (i+1)*nSym]);
        }
    }


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator) {}


  private:
    immutable bool[] _scMap;
    immutable size_t _nWLLearning;
    immutable size_t _nIter;
    immutable size_t _nCP, _nFFT, _nOS;

    FFTWObject!Complex _fftw;
    OverlapSaveRegenerator2!C _regenerator;

    immutable(C)[] _inputs;
    immutable(C)[] _desired;
    Complex!real _iqTX, _iqRX;
    Complex!real[P] _paCoefs;
    Complex!real[] _channelFreqResponse;
    size_t _nBufferedSymbols;

    C[] _buffer4Regen;


    void applyForEachSymbol(Flag!"learning" doLearning)(in C[] input, in C[] desired, C[] errors)
    {
        if(doLearning && _nBufferedSymbols < _nWLLearning) {
            _inputs ~= input;
            _desired ~= desired;
            ++_nBufferedSymbols;

            if(_nBufferedSymbols == _nWLLearning) {
                C[] ips = _inputs.dup;
                C[] dss = _desired.dup;

                C[][] freqX, freqY;
                estimateIQCoefs(ips, dss, freqX, freqY);
                toIQless(ips, dss, freqX, freqY);
                foreach(iIter; 0 .. _nIter){
                    estimateCFR(ips, dss);
                    estimatePACoefs(ips, dss);
                }
            }
        }

        _buffer4Regen[] = input[];
        foreach(ref e; _buffer4Regen) e = e + _iqTX * e.conj;
        foreach(ref e; _buffer4Regen){
            auto c = e;
            foreach(p; 1 .. P) e += _paCoefs[p] * c * c.sqAbs^^p;
        }
        auto state = RegeneratorMISOState(this, _nFFT * _nOS);
        C[][] ip2 = new C[][](input.length, 1);
        foreach(i; 0 .. input.length) ip2[i][0] = _buffer4Regen[i];
        auto es = errors;
        _regenerator.apply(state, ip2, es);
        assert(es is errors);

        foreach(i; 0 .. input.length)
            errors[i] = desired[i] - (errors[i] + errors[i].conj * _iqRX);
    }


    void fftAndRemoveHighOrder(in C[] input, in C[] desired, C[] freqY, bool isIQFree)
    {
        immutable nSym = (_nFFT + _nCP) * _nOS;

        _fftw.inputs!float[] = desired[];
        _fftw.fft!float();
        freqY[] = _fftw.outputs!float[];

        foreach(p; 1 .. P){
            _fftw.inputs!float[] = input[];
            if(!isIQFree) foreach(ref e; _fftw.inputs!float) e = e + _iqTX * e.conj;
            foreach(f, ref e; _fftw.inputs!float) e = _paCoefs[p] * e * (e.sqAbs^^p) * _channelFreqResponse[f];
            if(!isIQFree) foreach(ref e; _fftw.inputs!float) e = e + _iqRX * e.conj;
            _fftw.fft!float();
            auto ops = _fftw.outputs!float;
            foreach(f; 0 .. _nFFT * _nOS)
                freqY[f] -= ops[f];
        }
    }


    void estimateIQCoefs(in C[] input, in C[] desired, ref C[][] freqX_dst, ref C[][] freqY_dst)
    {
        immutable nSym = (_nFFT + _nCP) * _nOS;

        C[][] freqX;
        C[][] freqY;
        foreach(i; 0 .. _nWLLearning){
            _fftw.inputs!float[] = input[_nCP * _nOS + i * nSym .. (i+1) * nSym];
            _fftw.fft!float();
            freqX ~= _fftw.outputs!float.dup;

            freqY ~= new C[_nFFT * _nOS];
            fftAndRemoveHighOrder(input[_nCP * _nOS + i * nSym .. (i+1) * nSym], desired[_nCP * _nOS + i * nSym .. (i+1) * nSym], freqY[i], false);
        }
        freqX_dst = freqX;
        freqY_dst = freqY;

        // H_10(f), H_01(f)を推定する
        C[2][size_t] hlist;
        {
            Slice!(2, C*) mx = new C[2 * _nWLLearning].sliced(2, _nWLLearning);
            C[] y = new C[_nWLLearning];
            foreach(f; getSCIndex4IQ()){
                foreach(i; 0 .. _nWLLearning){
                    mx[0, i] = freqX[i][f];
                    mx[1, i] = freqX[i][f == 0 ? 0 : $-f].conj;
                    y[i] = freqY[i][f];
                }

                auto h = leastSquareEstimate(mx, y);
                hlist[f] = [h[0], h[1]];
            }
        }

        // H_10(f), H_01(f)から係数を推定する
        {
            Slice!(2, C*) mx = new C[2 * hlist.length].sliced(2, hlist.length);
            C[] y = new C[hlist.length];
            foreach(i, f; hlist.keys){
                assert(f != 0);
                immutable invFreq = _nFFT * _nOS - f;
                assert(invFreq in hlist);

                mx[0, i] = hlist[f][0];
                mx[1, i] = hlist[invFreq][0].conj;
                y[i] = hlist[f][1];
            }

            auto iqtxrx = leastSquareEstimate(mx, y);
            _iqTX = iqtxrx[0];
            _iqRX = iqtxrx[1];
        }
    }


    void toIQless(C[] input, C[] desired, C[][] freqX, C[][] freqY)
    {
        immutable nSym = (_nFFT + _nCP) * _nOS;

        foreach(i; 0 .. input.length) {
            input[i] = input[i] + _iqTX * input[i].conj;
            desired[i] = (desired[i] - _iqRX * desired[i].conj) / (1 - _iqRX.sqAbs);
        }

        foreach(i; 0 .. freqX.length) foreach(f; 0 .. _nFFT * _nOS / 2) {
            auto f1 = f;
            auto f2 = f == 0 ? 0 : (_nFFT * _nOS - f);
            auto x1 = freqX[i][f1] + _iqTX * freqX[i][f2].conj;
            auto x2 = freqX[i][f2] + _iqTX * freqX[i][f1].conj;
            auto y1 = (freqY[i][f1] - _iqRX * freqY[i][f2].conj) / (1 - _iqRX.sqAbs);
            auto y2 = (freqY[i][f2] - _iqRX * freqY[i][f1].conj) / (1 - _iqRX.sqAbs);
            freqX[i][f1] = x1;
            freqX[i][f2] = x2;
            freqY[i][f1] = y1;
            freqY[i][f2] = y2;
        }
    }


    void estimateCFR(in C[] input, in C[] desired)
    {
        immutable nSym = (_nFFT + _nCP) * _nOS;

        auto ips = _fftw.inputs!float;
        auto ops = _fftw.outputs!float;
        C[][] freqX = new C[][](_nWLLearning, _nFFT * _nOS);
        C[][] freqY = new C[][](_nWLLearning, _nFFT * _nOS);
        foreach(i; 0 .. _nWLLearning) foreach(f; 0 .. _nFFT * _nOS){
            freqX[i][f] = complexZero!C;
            freqY[i][f] = complexZero!C;
        }

        foreach(i; 0 .. _nWLLearning){
            _fftw.inputs!float[] = desired[_nCP * _nOS + i * nSym .. (i+1) * nSym];
            _fftw.fft!float();
            foreach(f; 0 .. _nFFT * _nOS) freqY[i][f] += ops[f];

            foreach(p; 0 .. P){
                _fftw.inputs!float[] = input[_nCP * _nOS + i * nSym .. (i+1) * nSym];
                foreach(ref e; ips) e = _paCoefs[p] * e * e.sqAbs^^p;
                _fftw.fft!float();
                foreach(f; 0 .. _nFFT * _nOS) freqX[i][f] += ops[f];
            }
        }

        foreach(f; 0 .. _nFFT * _nOS){
            Complex!real num = complexZero!C,
                         den = complexZero!C;

            foreach(i; 0 .. _nWLLearning){
                num += freqX[i][f].conj * freqY[i][f];
                den += freqX[i][f].sqAbs;
            }

            _channelFreqResponse[f] = num / den;
        }
    }


    void estimatePACoefs(in C[] input, in C[] desired)
    {
        immutable nLearningSymbols = 2;
        immutable nSym = (_nFFT + _nCP) * _nOS;

        // 前後二つのSWAPシンボルを加算する
        C[][][P] freqX;
        C[][] freqY;
        auto ips = _fftw.inputs!float;
        auto ops = _fftw.outputs!float;
        foreach(i; 0 .. nLearningSymbols){
            foreach(p; 0 .. P){
                _fftw.inputs!float[] = input[_nCP * _nOS + i * nSym .. (i+1) * nSym];
                foreach(ref e; ips) e = e * e.sqAbs^^p;
                _fftw.fft!float();
                if(i % 2 == 0)  freqX[p] ~= _fftw.outputs!float.dup;
                else foreach(f; 0 .. _nFFT * _nOS) freqX[p][i/2][f] += ops[f];
            }

            _fftw.inputs!float[] = desired[_nCP * _nOS + i * nSym .. (i+1) * nSym];
            _fftw.fft!float();
            if(i % 2 == 0)  freqY ~= _fftw.outputs!float.dup;
            else foreach(f; 0 .. _nFFT * _nOS) freqY[i/2][f] += ops[f];
        }

        // a_pを推定していく
        foreach_reverse(p; 1 .. P){
            Complex!real num = complexZero!C,
                         den = complexZero!C;

            foreach(f; getSCIndex4PthAMPCoef(p)) foreach(i; 0 .. nLearningSymbols / 2){
                auto x = freqX[p][i][f] * _channelFreqResponse[f];
                auto y = freqY[i][f];
                num += x.conj * y;
                den += x.sqAbs;
            }

            _paCoefs[p] = num / den;

            // 干渉除去する
            foreach(f; 0 .. _nFFT * _nOS) foreach(i; 0 .. nLearningSymbols / 2)
                freqY[i][f] -= freqX[p][i][f] * _channelFreqResponse[f] * _paCoefs[p];
        }
    }


    auto getSCIndex4IQ()
    {
        size_t n;
        foreach(i; 0 .. _nFFT * _nOS)
            if(_scMap[i]) ++n;

        n /= 2;
        n -= 10;
        auto r1 = iota(1, 1 + n);
        auto r2 = iota(_nFFT * _nOS - n, _nFFT * _nOS);
        return r1.chain(r2);
    }


    auto getSCIndex4PthAMPCoef(size_t p)
    {
        size_t n;
        foreach(i; 0 .. _nFFT * _nOS)
            if(_scMap[i]) ++n;

        n /= 2;
        auto r1 = iota(1+(p*2-1)*n, 1+(p*2+1)*n-10);
        auto r2 = iota(_nFFT * _nOS - (p*2+1)*n+10, _nFFT * _nOS - (p*2-1)*n);
        return r1.chain(r2);
    }


    static struct RegeneratorMISOState
    {
        this(IQInversionSuccessiveInterferenceCanceller c, size_t nFFT)
        {
            this.canceller = c;
            this.nFFT = nFFT;
        }

        IQInversionSuccessiveInterferenceCanceller canceller;
        size_t nFFT;

        void regenerate(C[][] distX, ref C[] dst)
        {
            if(dst.length != nFFT) dst.length = nFFT;

            foreach(f; 0 .. nFFT)
                dst[f] = canceller._channelFreqResponse[f] * distX[f][0];
        }
    }
}

