module dffdd.filter.sidelobe;

import std.algorithm;
import std.array;
import std.complex;
import std.json;
import std.range;
import std.typecons;

import mir.ndslice : slice, Slice, Contiguous;

import carbon.math : complexZero;

import dffdd.filter.freqdomain;
import dffdd.utils.fft;


final class SidelobeIterativeWLNL(C, size_t P)
{
    import mir.ndslice : Slice, Contiguous, sliced;
    import dffdd.utils.linalg;

    alias R = typeof(C.init.re);

    // 今のところP==2,　つまりx|x|^2の3次の項までしか考慮していない
    static assert(P == 2);


    this(size_t trainingSymbols, size_t nIter, size_t nFFT, size_t nCP, size_t nTone, size_t nOS)
    {
        _nTr = trainingSymbols;
        _nIter = nIter;
        _nSC = nTone;
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;

        _fftw = makeFFTWObject!Complex(_nFFT * _nOS);
        _regenerator = new OverlapSaveRegenerator2!C(1, _nFFT * _nOS);

        _iqTX = 0;
        _iqRX = 0;

        _paCoefs[0] = 1;
        _lnaCoefs[0] = 1;
        foreach(ref e; _paCoefs[1 .. $]) e = 0;
        foreach(ref e; _lnaCoefs[1 .. $]) e = 0;

        _channelFreqResponse = new Complex!double[_nFFT * _nOS];
        foreach(ref e; _channelFreqResponse) e = 0;

        _buffer4Regen = new C[(_nFFT + _nCP) * _nOS];
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


    JSONValue info()
    {
        static
        JSONValue cpxToJV(F)(Complex!F a){
            return JSONValue(["re": a.re, "im": a.im]);
        }

        JSONValue jv = [ "type": typeof(this).stringof ];
        jv["estTXIQCoef"] = cpxToJV(_iqTX);
        jv["estRXIQCoef"] = cpxToJV(_iqRX);
        jv["estPACoefs"] = _paCoefs[].map!cpxToJV().array();
        jv["estLNACoefs"] = _lnaCoefs[].map!cpxToJV().array();

        // 推定した周波数応答をインパルス応答に直す
        _fftw.inputs!double[] = _channelFreqResponse[];
        _fftw.ifft!double();
        jv["estImpResp"] = _fftw.outputs!double[].map!cpxToJV().array();

        return jv;
    }


  private:
    immutable size_t _nTr;
    immutable size_t _nIter;
    immutable size_t _nCP, _nFFT, _nOS, _nSC;

    FFTWObject!Complex _fftw;
    OverlapSaveRegenerator2!C _regenerator;

    immutable(C)[] _transmits;
    immutable(C)[] _receives;
    size_t _nBufferedSymbols;

    Complex!real _iqTX, _iqRX;
    Complex!real[P] _paCoefs;
    Complex!real[P] _lnaCoefs;
    Complex!double[] _channelFreqResponse;

    C[] _buffer4Regen;


    void applyForEachSymbol(Flag!"learning" doLearning)(in C[] transmits, in C[] receives, C[] outputs)
    {
        if(doLearning && _nBufferedSymbols < _nTr) {
            _transmits ~= transmits;
            _receives ~= receives;
            ++_nBufferedSymbols;

            // 学習用の信号が溜まったら，一度に学習する
            if(_nBufferedSymbols == _nTr) {
                foreach(iIter; 0 .. _nIter) {
                    estimateWLModel();
                    estimateNLCoefs();
                }
            }
        }

        // // 受信信号を出力にコピーして
        // // 受信IQミキサと受信アンプの歪みの逆特性を与える
        outputs[] = receives[];
        foreach(ref e; outputs) {
            e = (e - _iqRX * e.conj) / (1 - _iqRX.sqAbs);
            e -= _lnaCoefs[1] * e * e.sqAbs;
        }

        // 推定解を用いて送信信号に歪みを与えていく
        _buffer4Regen[] = transmits[];
        foreach(ref e; _buffer4Regen){
            e = e + _iqTX * e.conj;
            immutable c = e;
            foreach(p; 1 .. P) e += _paCoefs[p] * c * c.sqAbs^^p;
        }

        // 推定したチャネルに通す
        auto state = RegeneratorMISOState(this, _nFFT * _nOS);
        auto reconsted = new C[outputs.length];
        _regenerator.apply(state, _buffer4Regen.map!"[a]".array(), reconsted);

        // foreach(ref e; reconsted) {
        //     e += _lnaCoefs[1] * e * e.sqAbs;
        //     e += e.conj * _iqRX;
        // }

        // 干渉除去
        foreach(i, ref e; outputs)
            e -= reconsted[i];
    }


    C[] removeHighOrderNL() const
    {
        // IQミキサ
        void applyIQMixer(ref C[] s, Complex!real c) {
            foreach(ref e; s) e += c * e.conj;
        }

        // アンプ
        void applyAmplifierNL(ref C[] s, size_t p, Complex!real c) {
            foreach(ref e; s) e *= c * e.sqAbs^^p;
        }

        // チャネル (OWS法)
        void applyChannel(ref C[] s) {
            auto regen = new OverlapSaveRegenerator2!C(1, _nFFT * _nOS);
            auto state = RegeneratorMISOState(this, _nFFT * _nOS);
            regen.apply(state, s.map!"[a]".array(), s);
        }

        // 除去
        void applyRemove(C[] from, C[] src) {
            foreach(i, ref e; from) e -= src[i];
        }


        auto removed = _receives.dup;
        foreach(p; 1 .. P) {
            // Mixer -> TXAmp -> Channel -> Mixerの成分を消す
            auto signal = _transmits.dup;
            applyIQMixer(signal, _iqTX);
            applyAmplifierNL(signal, p, _paCoefs[p]);
            applyChannel(signal);
            applyIQMixer(signal, _iqRX);
            applyRemove(removed, signal);

            // Mixer -> Channel -> RXAmp -> Mixerの成分を消す
            signal = _transmits.dup;
            applyIQMixer(signal, _iqTX);
            applyChannel(signal);
            applyAmplifierNL(signal, p, _lnaCoefs[p]);
            applyIQMixer(signal, _iqRX);
            applyRemove(removed, signal);
        }

        return removed;
        // return _receives.dup;
    }


    void estimateWLModel()
    {
        auto removedNLY = removeHighOrderNL()[_nCP*_nOS-1 .. $];   // 先頭からnCP-1分は消す
        auto mx = slice!C(_nCP*_nOS*2, removedNLY.length);
        foreach(i; 0 .. removedNLY.length) {
            foreach(j; 0 .. _nCP*_nOS){
                // j=0でi+_nCP-1にするためi+_nCP-1-jとなる
                mx[j, i]      = _transmits[i+_nCP*_nOS-1-j];
                mx[j+_nCP*_nOS, i] = _transmits[i+_nCP*_nOS-1-j].conj;
            }
        }

        // 最小二乗法でWLのインパルス応答を求める
        // [0 .. _nCP]にxに対するh_1[k], [_nCP .. 2*_nCP]にconj(x)に対するh_2[k]
        auto estimatedWLH = leastSquareEstimate(mx, removedNLY);
        auto h1 = estimatedWLH[0 .. _nCP*_nOS];
        auto h2 = estimatedWLH[_nCP*_nOS .. _nCP*_nOS * 2];

        // h1とh2からIQImbalanceの係数_iqTX, _iqRXを推定
        auto iqTXRX = leastSquareEstimate2(h1, h1.map!conj, h2);
        _iqTX = iqTXRX[0];
        _iqRX = iqTXRX[1];

        // h1, h2, iqTX, iqRXからチャネルhを推定
        auto h = new Complex!double[_nCP*_nOS];
        foreach(i; 0 .. _nCP*_nOS)
            h[i] = (h1[i] - _iqTX.conj * _iqRX * h1[i].conj) / (1 - _iqTX.sqAbs * _iqRX.sqAbs);

        // インパルス応答hから周波数応答Hに変換
        _fftw.inputs!double[0 .. _nCP*_nOS] = h[0 .. $];
        _fftw.inputs!double[_nCP*_nOS .. $] = complexZero!(Complex!double);
        _fftw.fft!double();
        _channelFreqResponse[] = _fftw.outputs!double[];
    }


    void estimateNLCoefs()
    {
        immutable nSym = (_nFFT + _nCP) * _nOS;
        
        // OFDMシンボルを周波数領域に変換する
        // 初期化後，freqPsiXとfreqPsiYについてはp=0についてfreqPsiX[i][j][p]==NaN
        C[][][P] freqPsiX;
        C[][][P] freqPsiY;
        C[][] freqY;

        auto ips = _fftw.inputs!R;
        auto ops = _fftw.outputs!R;
        foreach(i; 0 .. _nTr) {
            foreach(p; 1 .. P) {
                // i番目のシンボルを抜き出して非線形変換
                ips[] = _transmits[_nCP * _nOS + i * nSym .. (i+1) * nSym];
                foreach(ref e; ips) e = e + e.conj * _iqTX;
                foreach(ref e; ips) e = e * e.sqAbs^^p;
                _fftw.fft!R();
                foreach(f, ref e; ops) e *= _channelFreqResponse[f];
                freqPsiX[p] ~= ops.dup;

                ips[] = _transmits[_nCP * _nOS + i * nSym .. (i+1) * nSym];
                foreach(ref e; ips) e = e + e.conj * _iqTX;
                _fftw.fft!R();
                foreach(f, ref e; ips) e = ops[f] * _channelFreqResponse[f];
                _fftw.ifft!R();
                foreach(k, ref e; ips) e = ops[k] * ops[k].sqAbs^^p;
                _fftw.fft!R();
                freqPsiY[p] ~= ops.dup;
            }

            ips[] = _receives[_nCP * _nOS + i * nSym .. (i+1) * nSym];
            foreach(ref e; ips) e = (e - _iqRX * e.conj) / (1 - _iqRX.sqAbs);
            _fftw.fft!R();
            freqY ~= ops.dup;
        }

        // 係数を推定していく
        foreach_reverse(p; 1 .. P) {
            Complex!R[] vecPsiX, vecPsiY, vecY;
            
            foreach(f; getSCIndex4PthAMPCoef(p)) foreach(i; 0 .. _nTr){
                vecPsiX ~= freqPsiX[p][i][f];
                vecPsiY ~= freqPsiY[p][i][f];
                vecY ~= freqY[i][f];
            }

            auto estCoefs = leastSquareEstimate2(vecPsiX, vecPsiY, vecY);

            _paCoefs[p] = estCoefs[0];
            _lnaCoefs[p] = estCoefs[1];

            // 干渉除去する
            foreach(f; 0 .. _nFFT * _nOS) foreach(i; 0 .. _nTr){
                freqY[i][f] -= freqPsiX[p][i][f] * _paCoefs[p];
                freqY[i][f] -= freqPsiY[p][i][f] * _lnaCoefs[p];
            }
        }
    }


    auto getSCIndex4PthAMPCoef(size_t p) const
    {
        size_t n = _nSC;
        n /= 2;
        auto r1 = iota(1+(p*2-1)*n, 1+(p*2+1)*n-10);
        auto r2 = iota(_nFFT * _nOS - (p*2+1)*n+10, _nFFT * _nOS - (p*2-1)*n);
        return r1.chain(r2);
    }


    static struct RegeneratorMISOState
    {
        this(const SidelobeIterativeWLNL c, size_t nFFT)
        {
            this.canceller = c;
            this.nFFT = nFFT;
        }

        const(SidelobeIterativeWLNL) canceller;
        size_t nFFT;

        void regenerate(C[][] distX, ref C[] dst)
        {
            if(dst.length != nFFT) dst.length = nFFT;

            foreach(f; 0 .. nFFT)
                dst[f] = canceller._channelFreqResponse[f] * distX[f][0];
        }
    }
}

unittest 
{
    import dffdd.mod.ofdm;
    import std.random;
    import std.stdio;
    import std.math : PI;

    // 1シンボル52トーン，FFTサイズ64, CPサイズ16のOFDM
    auto mod = new OFDM!(Complex!double)(64, 16, 52, 4);

    // 100シンボル分のBPSK
    auto bpsk = Random().map!(a => Complex!double(a%2 == 0 ? 1 : -1, 0)).take(52*100).array();
    Complex!double[] ofdm;
    mod.modulate(bpsk, ofdm);
    foreach(ref e; ofdm) e *= 10;

    Complex!double[] signal = ofdm.dup;

    immutable Complex!real trueTXIQImb = 0.05 * std.complex.expi(0.1*PI);
    foreach(ref e; signal)
        e += e.conj * trueTXIQImb;

    immutable Complex!real trueTXPACoef = 0.03 * std.complex.expi(0.4*PI);
    foreach(ref e; signal)
        e += e * e.sqAbs * trueTXPACoef;

    Random rnd;
    rnd.seed(114514);
    Complex!double[] channel = new Complex!double[64];
    foreach(i; 0 .. 64)
        channel[i] = std.complex.expi(uniform(-1.0f, 1.0f, rnd) * PI) / 8;

    {
        auto newsignal = signal.dup;
        newsignal[] = complexZero!(Complex!double);

        foreach(i; 64 .. signal.length)
            foreach(j; 0 .. 64)
            newsignal[i] += signal[i-j] * channel[j];
        
        signal = newsignal;
    }

    immutable Complex!real trueRXLNACoef = 0.03 * std.complex.expi(-0.1*PI);
    foreach(ref e; signal)
        e += e * e.sqAbs * trueRXLNACoef;

    immutable Complex!real trueRXIQImb = 0.05 * std.complex.expi(-0.4*PI);
    foreach(ref e; signal)
        e += e.conj * trueRXIQImb;

    alias Canceller = SidelobeIterativeWLNL!(Complex!double, 2);
    auto canceller = new Canceller(100, 1, 64, 16, 52, 4);
    canceller.apply!(Yes.learning)(ofdm, signal, signal.dup);

    real relErr(Complex!real x, Complex!real t)
    {
        return std.complex.abs(t - x) / std.complex.abs(t);
    }

    // // 相対誤差が0.1以下
    // writefln!"%s : %s"(canceller._iqTX, trueTXIQImb);
    // writefln!"%s : %s"(canceller._iqRX, trueRXIQImb);
    // writefln!"%s : %s"(canceller._paCoefs[1], trueTXPACoef);
    // writefln!"%s : %s"(canceller._lnaCoefs[1], trueRXLNACoef);

    canceller._fftw.inputs!double[] = canceller._channelFreqResponse[];
    canceller._fftw.ifft!double();
    // writeln(canceller._fftw.outputs!double[0 .. 8]);
    // writeln(channel[0 .. 8]);
    // assert(relErr(canceller._iqTX, trueTXIQImb) < 0.1);
    // assert(relErr(canceller._iqRX, trueRXIQImb) < 0.1);
    // assert(relErr(canceller._paCoefs[1], trueTXPACoef) < 0.1);
    // assert(relErr(canceller._lnaCoefs[1], trueRXLNACoef) < 0.1);

    import std.math : log10;

    auto dst = signal.dup;
    canceller.apply!(No.learning)(ofdm, signal, dst);
    dst = dst[$/4 .. $/4*3];
    auto rem = dst.fold!((a, b) => a + b.sqAbs)(0.0L);
    auto si = signal.fold!((a, b) => a + b.sqAbs)(0.0L);
    // writeln(10*log10(si/rem));
    // assert(10*log10(si/rem) > 70);  // 70 dB以上の除去量
}
