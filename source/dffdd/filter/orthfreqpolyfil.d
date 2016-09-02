module dffdd.filter.orthfreqpolyfil;


final class HammersteinFilterFor(alias BasisFuncs, alias genAdaptor)
{
    alias C = Complex!float;
    enum size_t P = BasisFuncs.length;

    this(in bool[] subcarieerMap, size_t taps, size_t nFFT, size_t nCP, size_t nOS)
    {
        _model = model;
        _fftw = makeFFTWObject(subcarieerMap.length);
        _taps = taps;
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;
        _isLastModeLearning = false;

        foreach(i, e; _scmap){
            if(e){
                _adaptor ~= genAdaptor(i);
                _state ~= new OneTapMultiFIRState!(Complex!float, BasisFuncs.length)();
            }else{
                _adaptor ~= null;
                _state ~= null;
            }
        }

        foreach(p; 0 .. P)
        {
            _lastTXSigs[p] = new Complex!float[subcarieerMap.length];
            _coefsOfFIRs[p] = new Complex!float[subcarieerMap.length];

            _lastTXSigs[] = Complex!float(0, 0);
            _coefsOfFIRs[] = Complex!float(0, 0);
        }

        adjustReserved(nFFT * nOS);
    }


    void apply(bool bLearning = true, C)(in C[] tx, in C[] rx, C[] output)
    in{
        assert(tx.length == rx.length);
        assert(tx.length == output.length);

      static if(bLearning)
      {
        assert(tx.length % (_OS * (_nFFT+_nCP)) == 0);

        // 先頭にCPがあるかチェック
        foreach(i; 0 .. _nOS * _nCP)
            assert(tx[i] == tx[i + _nOS * _nFFT]);
      }
    }
    body{
        if(bLearning)
        {
            immutable symLen = _nOS * (_nFFT * _nCP),
                      cpLen = _nOS * _nCP;

            foreach(sidx; 0 .. tx.length / symLen){
                const isymbol = tx[sidx * symLen .. (sidx+1)*symLen][cpLen .. $];
                const osymbol = rx[sidx * symLen .. (sidx+1)*symLen][cpLen .. $];
                C[][P] ffted;

                // 周波数領域での非線形送信信号の生成
                foreach(p, func; BasisFuncs){
                    ffted[p] = _reserved[p][0 .. _nOS * _nFFT];

                    foreach(i, ref e; _fftw.inputs!float[])
                        e = func(isymbol[i]);

                    _fftw.fft!float();
                    ffted[p][] = _fftw.outputs!float[];
                }

                // 受信信号の周波数成分
                _fftw.inputs!float[] = osymbol[];
                _fftw.fft!float();
                auto rxffted = _fftw.outputs!float;

                // 全周波数で，適応フィルタにかける
                foreach(freq; 0 .. nfft * ros) if(_state[freq] !is null) {
                    C[P] inps;
                    foreach(p, f; BasisFuncs)
                        inps[p] = ffted[p][freq];

                    _state[freq].update(inps);
                    C error = _state.error(rxffted[freq]);
                    _adaptor[freq].adapt(_state[freq], error);
                }
            }

            _isLastModeLearning = true;
        }
        else
        {
            if(_isLastModeLearning){
                // calculate coefs
                foreach(p, func; BasisFuncs){
                    auto ips = _fftw.inputs!float;

                    foreach(freq; 0 .. nfft * ros){
                        if(_state[freq] !is null)
                            ips[freq] = _state[freq].weight[0][p];
                        else
                            ips[freq] = 0;
                    }

                    _fftw.ifft!float();
                    auto ops = _fftw.outputs!float;
                    ops[_taps .. $] = Complex!float(0, 0);
                    assert(_taps <= ops.length / 2);

                    _fftw.inputs!float[] = _fftw.outputs!float[];
                    _fftw.fft!float();
                    _coefsOfFIRs[p][] = _fftw.outputs!float[];
                }

                _isLastModeLearning = false;
            }

            outputs[] = rx[];

            immutable batchLen = _fftw.inputs!float.length / 2;

            size_t totalLoopSize = tx.length / batchLen;
            if(tx.length % batchLen != 0)
                ++totalLoopSize;

            foreach(bidx; 0 .. totalLoopSize){
                immutable startIdx = batchLen * bidx,
                          remainLen = tx.length - startIdx,
                          thisBatchLen = min(remainLen, batchLen),
                          endIdx = startIdx + thisBatchLen;

                foreach(p, func; BasisFuncs){
                    auto ips = _fftw.inputs!float;
                    auto ops = _fftw.outputs!float;

                    immutable preLen = _lastTXSigs[p].length;

                    ips[0 .. preLen] = _lastTXSigs[p][];
                    ips[preLen .. preLen + thisBatchLen] = tx[startIdx .. endIdx];
                    ips[preLen + thisBatchLen .. $] = Complex!float(0, 0);
                    _fftw.fft!float();
                    ips[] = ops[];
                    foreach(i; 0 .. ips.length) ips[i] *= _coefsOfFIRs[p][i];
                    _fftw.ifft!float();
                    
                    foreach(i; 0 .. thisBatchLen)
                        output[startIdx + i] -= ops[preLen + i];

                    // update lastTXSigs
                    immutable copyLenFromNew = min(thisBatchLen, preLen);
                    immutable copyLenFromOld = copyLenFromNew < preLen ? (preLen - copyLenFromNew) : 0;
                    foreach(i; 0 .. preLen){
                        if(i < copyLenFromOld)
                            _lastTXSigs[i] = _lastTXSigs[$ - copyLenFromOld + i];
                        else
                            _lastTXSigs[i] = tx[endIdx - copyLenFromNew + (i - copyLenFromOld)];
                    }
                }
            }
        }
    }


  private:
    alias AdaptorType = typeof(genAdaptor(0));

    immutable size_t _taps, _nFFT, _nCP, _nOS;
    bool _isLastModeLearning;
    FFTWObject!(Complex!float) _fftw;
    OneTapMultiFIRState!(Complex!float, BasisFuncs.length)[] _state;
    AdaptorType[] _adaptor;
    Complex!float[][P] _coefsOfFIRs;
    Complex!float[][P] _lastTXSigs;

  static:
    C[][P] _reserved;

    static this()
    {

    }


    void adjustReserved(size_t size)
    {
        foreach(p; 0 .. P)
            if(_reserved[p].length < size)
                _reserved[p].length = size;
    }
}

unittest
{

}

