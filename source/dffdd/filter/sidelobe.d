module dffdd.filter.sidelobe;

import std.algorithm;
import std.array;
import std.complex;
import std.exception;
import std.json;
import std.range;
import std.typecons;

import mir.ndslice : slice, Slice, Contiguous;

import dffdd.filter.freqdomain;
import dffdd.utils.fft;


final class SidelobeIterativeWLNL(C, size_t P)
{
    import mir.ndslice : Slice, Contiguous, sliced;
    import dffdd.utils.linalg;

    alias R = typeof(C.init.re);

    // 今のところP==2,　つまりx|x|^2の3次の項までしか考慮していない
    // static assert(P == 2);


    this(size_t trainingSymbols, size_t nIter, size_t nFFT, size_t nCP, size_t nTone, size_t nOS, size_t nImpulseTaps,
        Flag!"isChFreqEst" isChFreqEst = Yes.isChFreqEst, Flag!"isInvertRX" isInvertRX = Yes.isInvertRX,
        Flag!"useNewton" useNewton = No.useNewton, size_t nNewtonIter = 10,
        Flag!"use3rdSidelobe" use3rdSidelobe = No.use3rdSidelobe,
        size_t nUseSCForEstNonlin = 50,
        string estOrder = "IHD")
    {
        _nTr = trainingSymbols;
        _nIter = nIter;
        _nNewtonIter = nNewtonIter;
        _nSC = nTone;
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;
        _nImpulseTaps = nImpulseTaps;
        _isChFreqEst = cast(bool)isChFreqEst;
        _isInvertRX = cast(bool)isInvertRX;
        _useNewton = cast(bool)useNewton;
        _use3rdSidelobe = cast(bool)use3rdSidelobe;
        _estOrder = estOrder;

        _nUseSC = min(nUseSCForEstNonlin, nTone);
        enforce(_nUseSC % 2 == 0);

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
    immutable size_t _nIter, _nNewtonIter, _nUseSC;
    immutable size_t _nCP, _nFFT, _nOS, _nSC;
    immutable size_t _nImpulseTaps;
    immutable bool _isChFreqEst, _isInvertRX, _useNewton, _use3rdSidelobe;
    immutable string _estOrder;

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
                void delegate()[char] stages = [
                    'I': &estimateIQImbalance,
                    'H': &estimateChannel,
                    'D': &estimateNLCoefs,
                ];

                foreach(iIter; 0 .. _nIter) {
                    foreach(char c; _estOrder)    // estOrderで指定された順に推定する
                        stages[c]();
                }
            }
        }

        // 受信信号を出力にコピーして
        outputs[] = receives[];

        // 受信IQミキサと受信アンプの歪みの逆特性を与える
        if(_isInvertRX){
            foreach(ref e; outputs) {
                e = (e - _iqRX * e.conj) / (1 - _iqRX.sqAbs);
                // e = invertPoly(e, C(_lnaCoefs[1]));
                e = invertPoly(e, _lnaCoefs[]);
            }
        }

        // 推定解を用いて送信信号に歪みを与えていく
        _buffer4Regen[] = transmits[];
        foreach(ref e; _buffer4Regen){
            e = e + _iqTX * e.conj;
            immutable c = e;
            // foreach(p; 1 .. P) e += _paCoefs[p] * c * c.sqAbs^^p;
            e = mypoly(e, _paCoefs[]);
        }

        // 推定したチャネルに通す
        auto state = RegeneratorMISOState(this, _nFFT * _nOS);
        auto reconsted = new C[outputs.length];
        _regenerator.apply(state, _buffer4Regen.map!"[a]".array(), reconsted);

        // 受信機側の非線形性を再現する
        if(!_isInvertRX){
            foreach(ref e; reconsted) {
                // e += _lnaCoefs[1] * e * e.sqAbs;
                e = mypoly(e, _lnaCoefs[]);
                e += e.conj * _iqRX;
            }
        }

        // 干渉除去
        foreach(i, ref e; outputs)
            e -= reconsted[i];
    }


    C[][2] removeHighOrderNL() const
    {
        auto txs = _transmits.dup;
        foreach(ref e; txs) {
            e += e.conj * _iqTX;
            // e += _paCoefs[1] * e * e.sqAbs;
            e = mypoly(e, _paCoefs[]);
            e = (e - _iqTX * e.conj) / (1 - _iqTX.sqAbs);
        }

        auto rxs = _receives.dup;
        foreach(ref e; rxs) {
            e = (e - _iqRX * e.conj) / (1 - _iqRX.sqAbs);
            e = invertPoly(e, _lnaCoefs[]);
            e += e.conj * _iqRX;
        }

        return [txs, rxs];
        // return _receives.dup;
    }


    void estimateChannel()
    {
        import std.stdio;
        immutable nSym = (_nFFT + _nCP) * _nOS;
        auto removedNLXY = removeHighOrderNL();
        auto txs = removedNLXY[0];
        auto rxs = removedNLXY[1];        // 先頭からCPは消す

        foreach(ref e; txs)
            e += e.conj * _iqTX;
        
        foreach(ref e; rxs)
            e = (e - e.conj * _iqRX) / (1 -_iqRX.sqAbs);

        if(_isChFreqEst) {
            txs = txs[_nCP*_nOS .. $];
            rxs = rxs[_nCP*_nOS .. $];

            C[][] txsFreq;
            C[][] rxsFreq;
            foreach(i; 0 .. _nTr) {
                _fftw.inputs!R[] = txs[i*nSym .. i*nSym + _nFFT*_nOS];
                _fftw.fft!R();
                txsFreq ~= _fftw.outputs!R.dup;

                _fftw.inputs!R[] = rxs[i*nSym .. i*nSym + _nFFT*_nOS];
                _fftw.fft!R();
                rxsFreq ~= _fftw.outputs!R.dup;
            }

            auto freqs = (){
                auto fpos = iota(1, _nSC/2+1);
                auto fs = fpos.array();
                foreach(f; fpos)
                    fs ~= _nFFT*_nOS - f;
                return fs;
            }();

            // import std.stdio;
            // writeln(freqs);

            C[] estFreqResp;
            foreach(k; freqs) {
                C corrSum = C(0, 0);
                C sqSum = C(0, 0);
                foreach(i; 0 .. _nTr) {
                    corrSum += rxsFreq[i][k] * txsFreq[i][k].conj;
                    sqSum += txsFreq[i][k].sqAbs;
                }

                estFreqResp ~= corrSum / sqSum;
            }

            auto impResp = estimateImpulseResponseFromFrequencyResponse(_nFFT*_nOS, _nImpulseTaps, estFreqResp, freqs, 0);
            foreach(i, ref e; _fftw.inputs!double[0 .. _nImpulseTaps])
                e = impResp[i];
            _fftw.inputs!double[_nImpulseTaps .. $] = Complex!double(0);
            _fftw.fft!double();
            _channelFreqResponse[] = _fftw.outputs!double[];
        }else{
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

            _fftw.inputs!double[_nImpulseTaps .. $] = Complex!double(0);
            _fftw.fft!double();
            _channelFreqResponse[] = _fftw.outputs!double[];
        }
    }


    void estimateIQImbalance()
    {
        immutable nSym = (_nFFT + _nCP) * _nOS;
        auto removedNLXY = removeHighOrderNL();
        auto txs = removedNLXY[0][_nCP*_nOS .. $];
        auto rxs = removedNLXY[1][_nCP*_nOS .. $];        // 先頭からCPは消す


        C[][] freqTX;
        C[][] freqRX;
        foreach(i; 0 .. _nTr) {
            _fftw.inputs!R[] = txs[i*nSym .. i*nSym + _nFFT * _nOS];
            _fftw.fft!R();
            freqTX ~= _fftw.outputs!R.dup;

            _fftw.inputs!R[] = rxs[i*nSym .. i*nSym + _nFFT * _nOS];
            _fftw.fft!R();
            freqRX ~= _fftw.outputs!R.dup;
        }

        // H_10(f), H_01(f)を推定する
        C[2][size_t] hlist;
        {
            foreach(f; getSCIndex4IQ()){
                immutable invf = (f == 0 ? 0 : (_nFFT * _nOS)-f);

                auto h = leastSquareEstimate2(freqTX.map!(a => a[f]),
                                              freqTX.map!(a => a[invf].conj),
                                              freqRX.map!(a => a[f]));
                
                hlist[f] = h;
            }
        }

        // H_10(f), H_01(f)から係数を推定する
        {
            auto ks = hlist.keys;
            auto invks = ks.zip(repeat(_nFFT * _nOS)).map!"a[0] == 0 ? 0 : a[1] - a[0]";
            auto iqtxrx = leastSquareEstimate2(ks.zip(hlist.repeat).map!"a[1][a[0]][0]",
                                               invks.zip(hlist.repeat).map!(a => a[1][a[0]][0].conj),
                                               ks.zip(hlist.repeat).map!"a[1][a[0]][1]");
            _iqTX = iqtxrx[0];
            _iqRX = iqtxrx[1];
        }
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
            foreach(p; 0 .. P) {
                // i番目のシンボルを抜き出して非線形変換
                ips[] = _transmits[_nCP * _nOS + i * nSym .. (i+1) * nSym];
                foreach(ref e; ips) e = e + e.conj * _iqTX;
                foreach(ref e; ips) e = e * e.sqAbs^^p;
                _fftw.fft!R();
                foreach(f, ref e; ops) e *= _channelFreqResponse[f];
                freqPsiX[p] ~= ops.dup;

                ips[] = _transmits[_nCP * _nOS + i * nSym .. (i+1) * nSym];
                foreach(ref e; ips) e = e + e.conj * _iqTX;     // apply TX I/Q imbalance
                foreach(ref e; ips) e = mypoly(e, _paCoefs[]);  // apply TX PA Distortion
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
            // foreach(f, ref e; ops) e -= freqPsiX[0][i][f];      // 線形成分を抜いておく
            freqY ~= ops.dup;
        }

        // メインローブを使うなら
        if(!_use3rdSidelobe) {
            // 線形成分を抜いておく
            foreach(i; 0 .. _nTr)
                foreach(f, ref e; freqY[i])
                    e -= freqPsiX[0][i][f]; 
        }

        size_t[] targetFS;
        if(_use3rdSidelobe)
            targetFS = getSCIndex4PthAMPCoef(1).array;  // 3次IMDを利用
        else
            targetFS = getSCIndex4IQ(Yes.forNL).array;  // メインローブを利用

        // if(_use3rdSidelobe)
        {
            import mir.ndslice;
            import kaleidic.lubeck;

            Complex!R[] vecY = new Complex!R[](targetFS.length * _nTr);
            auto matX = slice!(Complex!R)(targetFS.length * _nTr, (P-1) * 2);

            foreach(ri, idxFI; cartesianProduct(targetFS, iota(_nTr)).enumerate) {
                size_t f = idxFI[0];        // 周波数インデックス
                size_t i = idxFI[1];        // シンボルインデックス

                vecY[ri] = freqY[i][f];

                foreach(p; 1 .. P){
                    matX[ri, (p-1)]         = freqPsiX[p][i][f];
                    matX[ri, (p-1) + (P-1)] = freqPsiY[p][i][f];
                }
            }

            auto coefs = matX.transposed.slice.leastSquareEstimateColumnMajor(vecY);
            foreach(p; 1 .. P) {
                _paCoefs[p] = coefs[(p-1) ];
                _lnaCoefs[p] = coefs[(p-1) + (P-1)];
            }
        }
        version(none)
        {
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
    }


    auto getSCIndex4IQ(Flag!"forNL" forNL = No.forNL)
    {
        size_t n = _nSC;

        n /= 2;
        auto r1 = iota(1, 1 + n);
        auto r2 = iota(_nFFT * _nOS - n, _nFFT * _nOS).retro();

        if(forNL)
            return r1[0 .. _nUseSC/2].chain(r2[0 .. _nUseSC/2]);
        else
            return r1.chain(r2);
    }


    auto getSCIndex4PthAMPCoef(size_t p) const
    {
        size_t n = _nSC;
        n /= 2;
        auto r1 = iota(1+(p*2-1)*n, 1+(p*2+1)*n)[0 .. _nUseSC/2];
        auto r2 = iota(_nFFT * _nOS - (p*2+1)*n, _nFFT * _nOS - (p*2-1)*n).retro()[0 .. _nUseSC / 2];
        return r1.chain(r2);
    }


    static
    C mypoly(C2)(C x, in C2[] coefs)
    {
        C sum = C(0);
        C xp = x;
        real sq = x.sqAbs;
        foreach(i, e; coefs) {
            sum += xp * e;
            xp *= sq;
        }

        return sum;
    }


    C invertPoly(C2)(C y, in C2[] coefs) const
    {
        if(_useNewton) {
            real yabs = y.abs;
            real e = yabs;
            foreach(i; 0 .. _nNewtonIter) {
                if(e == 0)
                    return C(0, 0);

                real fx = mypoly(C(e), coefs).abs - yabs;

                real sumR = 0, sumI = 0, sumRD = 0, sumID = 0;
                foreach(p, c; coefs) {
                    immutable ep = e^^(2*p);
                    sumR += c.re * ep * e;
                    sumI += c.im * ep * e;
                    sumRD += (2*p+1) * c.re * ep;
                    sumID += (2*p+1) * c.im * ep;
                }

                import std.math : sqrt;
                real dfx = 1 / sqrt(sumR^^2 + sumI^^2);
                dfx *= sumR * sumRD + sumI * sumID;

                e -= fx / dfx;
            }

            immutable real xabs = e;
            immutable C ypoly = mypoly(C(xabs), coefs);
            immutable C ytheta = ypoly / ypoly.abs;
            return C(xabs * (y/yabs) / ytheta);
        } else {
            return C(y - coefs[1] * y * y.sqAbs);
        }
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

// unittest 
// {
//     import dffdd.mod.ofdm;
//     import std.random;
//     import std.stdio;
//     import std.math : PI;

//     // 1シンボル52トーン，FFTサイズ64, CPサイズ16のOFDM
//     auto mod = new OFDM!(Complex!double)(64, 16, 52, 4);

//     // 100シンボル分のBPSK
//     auto bpsk = Random().map!(a => Complex!double(a%2 == 0 ? 1 : -1, 0)).take(52*100).array();
//     Complex!double[] ofdm;
//     mod.modulate(bpsk, ofdm);
//     foreach(ref e; ofdm) e *= 10;

//     Complex!double[] signal = ofdm.dup;

//     immutable Complex!real trueTXIQImb = 0.05 * std.complex.expi(0.1*PI);
//     foreach(ref e; signal)
//         e += e.conj * trueTXIQImb;

//     immutable Complex!real trueTXPACoef = 0.03 * std.complex.expi(0.4*PI);
//     foreach(ref e; signal)
//         e += e * e.sqAbs * trueTXPACoef;

//     Random rnd;
//     rnd.seed(114514);
//     Complex!double[] channel = new Complex!double[64];
//     foreach(i; 0 .. 64)
//         channel[i] = std.complex.expi(uniform(-1.0f, 1.0f, rnd) * PI) / 8;

//     {
//         auto newsignal = signal.dup;
//         newsignal[] = Complex!double(0);

//         foreach(i; 64 .. signal.length)
//             foreach(j; 0 .. 64)
//             newsignal[i] += signal[i-j] * channel[j];
        
//         signal = newsignal;
//     }

//     immutable Complex!real trueRXLNACoef = 0.03 * std.complex.expi(-0.1*PI);
//     foreach(ref e; signal)
//         e += e * e.sqAbs * trueRXLNACoef;

//     immutable Complex!real trueRXIQImb = 0.05 * std.complex.expi(-0.4*PI);
//     foreach(ref e; signal)
//         e += e.conj * trueRXIQImb;

//     alias Canceller = SidelobeIterativeWLNL!(Complex!double, 2);
//     auto canceller = new Canceller(100, 1, 64, 16, 52, 4);
//     canceller.apply!(Yes.learning)(ofdm, signal, signal.dup);

//     real relErr(Complex!real x, Complex!real t)
//     {
//         return std.complex.abs(t - x) / std.complex.abs(t);
//     }

//     // // 相対誤差が0.1以下
//     // writefln!"%s : %s"(canceller._iqTX, trueTXIQImb);
//     // writefln!"%s : %s"(canceller._iqRX, trueRXIQImb);
//     // writefln!"%s : %s"(canceller._paCoefs[1], trueTXPACoef);
//     // writefln!"%s : %s"(canceller._lnaCoefs[1], trueRXLNACoef);

//     canceller._fftw.inputs!double[] = canceller._channelFreqResponse[];
//     canceller._fftw.ifft!double();
//     // writeln(canceller._fftw.outputs!double[0 .. 8]);
//     // writeln(channel[0 .. 8]);
//     // assert(relErr(canceller._iqTX, trueTXIQImb) < 0.1);
//     // assert(relErr(canceller._iqRX, trueRXIQImb) < 0.1);
//     // assert(relErr(canceller._paCoefs[1], trueTXPACoef) < 0.1);
//     // assert(relErr(canceller._lnaCoefs[1], trueRXLNACoef) < 0.1);

//     import std.math : log10;

//     auto dst = signal.dup;
//     canceller.apply!(No.learning)(ofdm, signal, dst);
//     dst = dst[$/4 .. $/4*3];
//     auto rem = dst.fold!((a, b) => a + b.sqAbs)(0.0L);
//     auto si = signal.fold!((a, b) => a + b.sqAbs)(0.0L);
//     // writeln(10*log10(si/rem));
//     // assert(10*log10(si/rem) > 70);  // 70 dB以上の除去量
// }
