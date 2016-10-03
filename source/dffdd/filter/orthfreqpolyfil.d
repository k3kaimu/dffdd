module dffdd.filter.orthfreqpolyfil;

import std.algorithm;
import std.complex;
import std.stdio;

import dffdd.utils.fft;
import dffdd.filter.state;

final class FrequencyHammersteinFilter(alias genAdaptor, BasisFuncs...)
{
    alias C = Complex!float;
    enum size_t P = BasisFuncs.length;

    this(in bool[] subcarieerMap, size_t taps, size_t nFFT, size_t nCP, size_t nOS)
    in{
        //assert(subcarieerMap.length == nFFT * nOS);
    }
    body{
        //_model = model;
        _fftw = makeFFTWObject(nFFT * nOS);
        //_taps = taps;
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;
        _scMap = subcarieerMap.dup;

        foreach(i; 0 .. nFFT * nOS){
            _state ~= new OneTapMultiFIRFilterState!(Complex!float, BasisFuncs.length)();
            _adaptor ~= genAdaptor(i, subcarieerMap[i], _state[$-1]);
        }

        foreach(p; 0 .. P)
        {
            _lastTXSigs[p] = new Complex!float[nFFT * nOS/2];
            _coefsOfFIRs[p] = new Complex!float[nFFT * nOS];

            _lastTXSigs[p][] = Complex!float(0, 0);
            _coefsOfFIRs[p][] = Complex!float(0, 0);
        }

        adjustReserved(nFFT * nOS);
    }


    void apply(bool bLearning = true, C)(in C[] tx, in C[] rx, C[] output)
    in{
        assert(tx.length == rx.length);
        assert(tx.length == output.length);

      static if(bLearning)
      {
        assert(tx.length % (_nOS * (_nFFT+_nCP)) == 0);

        // 先頭にCPがあるかチェック
        foreach(i; 0 .. _nOS * _nCP)
            assert(tx[i] == tx[i + _nOS * _nFFT]);
      }
    }
    body{
        immutable symLen = _nOS * (_nFFT + _nCP),
                  cpLen = _nOS * _nCP;

        if(tx.length / symLen != 1){
            foreach(batchIndex; 0 .. tx.length / symLen)
                apply!bLearning(tx[symLen * batchIndex .. symLen * (batchIndex + 1)],
                                rx[symLen * batchIndex .. symLen * (batchIndex + 1)],
                                output[symLen * batchIndex .. symLen * (batchIndex + 1)]);

            return;
        }


        if(bLearning)
        {
            //immutable symLen = _nOS * (_nFFT + _nCP),
            //          cpLen = _nOS * _nCP;

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
                    //writeln(p, " : ", ffted[p]);
                }

                // 受信信号の周波数成分
                _fftw.inputs!float[] = osymbol[];
                _fftw.fft!float();
                auto rxffted = _fftw.outputs!float;

                // 全周波数で，適応フィルタにかける
                foreach(freq; 0 .. _nFFT * _nOS) {
                    //if(!isSwapped && !_scMap[freq]) continue;
                    C[P] inps;

                    foreach(p, f; BasisFuncs)
                        inps[p] = ffted[p][freq];

                    _state[freq].update(inps);
                    C error = _state[freq].error(rxffted[freq]);
                    _adaptor[freq].adapt(_state[freq], error);
                }
            }

            output[] = rx[];
        }
        //else
        {
            output[] = rx[];

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
                    //import std.stdio;
                    //writeln(preLen, " : ", thisBatchLen, " : ", ips.length);
                    foreach(i; 0 .. thisBatchLen)
                        ips[preLen + i] = func(tx[startIdx + i]);

                    ips[preLen + thisBatchLen .. $] = Complex!float(0, 0);
                    _fftw.fft!float();
                    ips[] = ops[];
                    foreach(i; 0 .. ips.length) ips[i] *= _state[i].weight[0][p];
                    _fftw.ifft!float();

                    foreach(i; 0 .. thisBatchLen)
                        output[startIdx + i] -= ops[preLen + i];

                    // update lastTXSigs
                    immutable copyLenFromNew = min(thisBatchLen, preLen);
                    immutable copyLenFromOld = copyLenFromNew < preLen ? (preLen - copyLenFromNew) : 0;
                    foreach(i; 0 .. preLen){
                        if(i < copyLenFromOld)
                            _lastTXSigs[p][i] = _lastTXSigs[p][$ - copyLenFromOld + i];
                        else
                            _lastTXSigs[p][i] = func(tx[endIdx - copyLenFromNew + (i - copyLenFromOld)]);
                    }
                }
            }
        }
    }


  private:
    alias AdaptorType = typeof(genAdaptor(0, true, new OneTapMultiFIRFilterState!(Complex!float, BasisFuncs.length)()));

    immutable size_t /*_taps,*/ _nFFT, _nCP, _nOS;
    immutable(bool[]) _scMap;
    FFTWObject!(Complex) _fftw;
    OneTapMultiFIRFilterState!(Complex!float, BasisFuncs.length)[] _state;
    AdaptorType[] _adaptor;
    Complex!float[][P] _coefsOfFIRs;
    Complex!float[][P] _lastTXSigs;

  static:
    C[][P] _reserved;

    void adjustReserved(size_t size)
    {
        foreach(p; 0 .. P)
            if(_reserved[p].length < size)
                _reserved[p].length = size;
    }
}
