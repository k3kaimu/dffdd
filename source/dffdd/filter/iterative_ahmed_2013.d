module dffdd.filter.iterative_ahmed_2013;



final class IterativeAhmed2013(C)
{
    this(size_t trainingSymbols, size_t nIter, size_t nIterNL, size_t nFFT, size_t nCP, size_t nTone, size_t nOS, size_t nImpulseTaps)
    {
        _nTr = trainingSymbols;
        _nIter = nIter;
        _nIterNL = nIterNL;
        _nNewtonIter = nNewtonIter;
        _nSC = nTone;
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;
        _nImpulseTaps = nImpulseTaps;
        _nBufferedSymbols = 0;

        _channelFreqResponse = new Complex!double[_nFFT * _nOS];
        foreach(ref e; _channelFreqResponse) e = 0;

        foreach(ref es; _buffer4Regen)
            es = new C[(_nFFT + _nCP) * _nOS];

        _fftw = makeFFTWObject!Complex(_nFFT * _nOS);
        _regenerator = new OverlapSaveRegenerator2!C(1, _nFFT * _nOS);
    }


    size_t inputBlockLength() @property
    {
        return (_nFFT + _nCP) * _nOS;
    }


    void apply(Flag!"learning" doLearning)(in C[] transmits, in C[] receives, C[] outputs)
    {
        immutable nSym = (_nFFT + _nCP) * _nOS;

        foreach(i; 0 .. transmits.length / nSym)
        {
            applyForEachSymbol!doLearning(transmits[i*nSym .. (i+1)*nSym], receives[i*nSym .. (i+1)*nSym], outputs[i*nSym .. (i+1)*nSym]);
        }
    }


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator) {}


  private:
    immutable size_t _nTr;
    immutable size_t _nIter, _nIterNL;
    immutable size_t _nCP, _nFFT, _nOS, _nSC;
    immutable size_t _nImpulseTaps;
    size_t _nBufferedSymbols;

    C _paCoef, _lnaCoef;
    C[] _channelFreqResponse;
    C[] _transmits, _receives;
    C[][2] _buffer4Regen;

    FFTWObject!Complex _fftw;
    OverlapSaveRegenerator2!C _regenerator;



    void applyForEachSymbol(Flag!"doLearning" doLearning)(in C[] transmits, in C[] receives, C[] outputs)
    {
        if(doLearning && _nBufferedSymbols < _nTr) {
            _transmits ~= transmits;
            _receives ~= receives;
            ++_nBufferedSymbols;

            // 学習用の信号が溜まったら，一度に学習する
            if(_nBufferedSymbols == _nTr) {
                // チャネルと非線形係数を交互に推定する
                foreach(iIter; 0 .. _nIter) {
                    estimateChannel();
                    estimateNLCoefs();
                }
            }
        }

        // 受信信号を出力にコピーして
        outputs[] = receives[];

        // 推定解を用いて送信信号にPAの歪みを与えていく
        _buffer4Regen[0] = transmits[];
        foreach(ref e; _buffer4Regen[0]){
            immutable c = e;
            e += _paCoefs * c * c.sqAbs;
        }

        // チャネルに通す
        auto state = RegeneratorMISOState(this, _nFFT * _nOS);
        _regenerator.apply(state, _buffer4Regen.map!"[a]".array(), _buffer4Regen[1]);

        // LNAの歪みを加える
        foreach(ref e; _buffer4Regen[1]) {
            immutable c = e;
            e += _lnaCoef * c * c.sqAbs;
        }

        // 干渉除去
        foreach(i, ref e; outputs)
            e -= _buffer4Regen[1][i];
    }


    void estimateChannel()
    {
        immutable nSym = (_nFFT + _nCP) * _nOS;
        auto txs = _transmits.dup;
        auto rxs = _receives.dup;

        rxs = rxs[_nImpulseTaps-1 .. $];
        auto mx = slice!C(_nImpulseTaps, rxs.length);
        foreach(i; 0 .. rxs.length) {
            foreach(j; 0 .. _nImpulseTaps){
                mx[j, i] = txs[i+j];
            }
        }

        auto estimatedH = leastSquareEstimateColumnMajor(mx, rxs);
        estimatedH.reverse();
        foreach(i, ref e; _fftw.inputs!double[0 .. _nImpulseTaps])
            e = estimatedH[i];

        _fftw.inputs!double[_nImpulseTaps .. $] = complexZero!(Complex!double);
        _fftw.fft!double();
        _channelFreqResponse[] = _fftw.outputs!double[];
    }


    void estimateNLCoefs()
    {
        
    }
}
