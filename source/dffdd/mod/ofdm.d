module dffdd.mod.ofdm;

import carbon.math;

import dffdd.utils.fft;

final class OFDM(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;
    alias RealType = typeof(C.init.re);

    this(uint nFFT, uint nCp, uint nTone, uint nUpSampling = 1)
    {
        _nFFT = nFFT;
        _nCp = nCp;
        _nTone = nTone;
        _nUpSampling = nUpSampling;
        _inpBuffer = new C[nFFT * nUpSampling];
        _fftw = makeFFTWObject!C(nFFT * nUpSampling);
    }


    size_t symInputLength() const @property { return _nTone; }
    size_t symOutputLength() const @property { return (_nFFT + _nCp) * _nUpSampling; }


    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    in{
        assert(inputs.length % this.symInputLength == 0);
    }
    body{
        _inpBuffer[] = complexZero!C;

        if(outputs.length != inputs.length / this.symInputLength * this.symOutputLength)
            outputs.length = inputs.length / this.symInputLength * this.symOutputLength;

        auto mainLobeH = _inpBuffer[0 .. _nFFT/2];
        auto mainLobeL = _inpBuffer[$ - _nFFT/2 .. $];

        foreach(i; 0 .. inputs.length / _nTone){
            //_inpBuffer[] = complexZero!C;
            mainLobeL[] = complexZero!C;
            mainLobeH[] = complexZero!C;

            auto inpTones = inputs[i*_nTone .. (i+1)*_nTone];
            if(_nTone == _nFFT){
                mainLobeH[0 .. _nTone/2] = inpTones[0 .. _nTone/2];
                mainLobeL[$-_nTone/2 .. $] = inpTones[_nTone/2 .. $];

            }else if(_nTone == _nFFT -1){
                mainLobeH[1 .. $] = inpTones[0 .. _nFFT/2-1];
                mainLobeL[0 .. $] = inpTones[_nFFT/2-1 .. _nTone];
            }else if(_nTone % 2 == 0){
                mainLobeH[1 .. _nTone/2 + 1] = inpTones[0 .. _nTone/2];
                mainLobeL[$ - _nTone/2 .. $] = inpTones[_nTone/2 .. _nTone];
            }else{
                mainLobeH[1 .. (_nTone+1)/2+1] = inputs[0 .. (_nTone+1)/2];
                mainLobeL[$ - _nTone/2 .. $] = inpTones[(_nTone+1)/2 .. _nTone];
            }

            auto dst = outputs[i*symOutputLength .. (i+1)*symOutputLength];
            .ifft!RealType(_fftw, _inpBuffer, dst[_nUpSampling * _nCp .. _nUpSampling * (_nCp + _nFFT)]);
            dst[0 .. _nUpSampling * _nCp] = dst[$ - _nUpSampling * _nCp .. $];

          version(DFFDD_OFDM_SIDELOBE_SUPPRESS)
          {
            foreach(j, ref e; dst[$-_nUpSampling .. $])
                e /= (j+2);
            foreach(j, ref e; dst[0 .. _nUpSampling])
                e /= (_nUpSampling+1 - j);
          }
            //dst[0] /= 2;
            //dst[1] /= 2;
            //dst[$-2] /= 2;
            //dst[$-1] /= 2;
            assert(dst.ptr == outputs.ptr + i*symOutputLength);
        }

        return outputs;
    }


    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    in{
        assert(inputs.length % this.symOutputLength == 0);
    }
    body{
        //outputs[] = complexZero!C;

        auto mainLobeH = _inpBuffer[0 .. _nFFT/2];
        auto mainLobeL = _inpBuffer[$ - _nFFT/2 .. $];

        if(outputs.length != inputs.length / this.symOutputLength * this.symInputLength)
            outputs.length = inputs.length / this.symOutputLength * this.symInputLength;

        foreach(i; 0 .. inputs.length / this.symOutputLength){
            auto inputSymbol = inputs[i*this.symOutputLength .. (i+1)*this.symOutputLength];
            auto outputSymbol = outputs[i*this.symInputLength .. (i+1)*this.symInputLength];
            .fft!RealType(_fftw, inputSymbol[_nUpSampling * _nCp .. $], _inpBuffer);

            if(_nTone == _nFFT){
                outputSymbol[] = _inpBuffer[];
            }else if(_nTone == _nFFT -1){
                outputSymbol[0 .. _nFFT/2 -1] = mainLobeH[1 .. $];
                outputSymbol[_nFFT/2 -1.. _nTone] = mainLobeL[0 .. $];
            }else if(_nTone % 2 == 0){
                outputSymbol[0 .. _nTone/2] = mainLobeH[1 .. _nTone/2 +1];
                outputSymbol[_nTone/2 .. _nTone] = mainLobeL[$ - _nTone/2 .. $];
            }else{
                outputSymbol[0 .. (_nTone+1)/2] = mainLobeH[1 .. (_nTone+1)/2+1];
                outputSymbol[(_nTone+1)/2 .. _nTone] = mainLobeH[$ - _nTone/2 .. $];
            }
        }

        return outputs;
    }


  private:
    FFTWObject!C _fftw;
    C[] _inpBuffer;
    uint _nFFT;
    uint _nCp;
    uint _nTone;
    uint _nUpSampling;
}

//
unittest
{
    import std.complex, std.math;

    // FFTサイズ: 8
    // サイクリックプレフィックス: 3
    // 使用サブキャリア数: 4
    // アップサンプリング率: 2
    auto ofdmMod = new OFDM!(Complex!float)(8, 3, 4, 2);

    Complex!float[] inps = [Complex!float(1, 1), Complex!float(-1, -1), Complex!float(1, 0), Complex!float(0, 1)];
    Complex!float[] res;

    // 変調
    ofdmMod.modulate(inps, res);

    assert(res.length == (8 + 3) * 2);

    // CPのチェック
    foreach(i; 0 .. 3 * 2)
        assert(res[i] == res[$ - 3*2 + i]);

    ofdmMod.demodulate(res.dup, res);

    assert(res.length == inps.length);

    foreach(i; 0 .. inps.length){
        assert(approxEqual(inps[i].re, res[i].re));
        assert(approxEqual(inps[i].im, res[i].im));
    }
}


/**
OFDM信号をオーバーサンプリングすることで折り返し雑音を再現します．
結果は周波数領域です．

Q: オーバーサンプリング率
BasisFuncs: 基底関数のリスト
tx: 送信レプリカ信号, N個のシンボルを一度に変換したい場合，tx.length == N * (numOfFFT + numOfCp)
numOfFFT: OFDM信号の実効FFTサイズ
numOfCp : OFDM信号の実効CPサイズ

Return: is(typeof(result[t][p][f]) == typeof(tx[0]))であり，
        result[t][p][f]は，tシンボル目のpバンド目の周波数f成分となる
*/
template generateOFDMAliasSignalAtFrequency(size_t Q, BasisFuncs...)
{
    import std.algorithm;
    import std.math;

    import carbon.complex;
    import dffdd.utils.fft;

    enum size_t QX = Q % 2 == 0 ? Q+1 : Q;

    C[][BasisFuncs.length * QX][] generateOFDMAliasSignalAtFrequency(C)(in C[] tx, size_t numOfFFT, size_t numOfCp)
    if(isComplex!C)
    in{
        assert(tx.length % (numOfFFT + numOfCp) == 0);
    }
    body{
        alias Cpx = complexTypeTemplate!C;
        alias F = typeof(C.init.re);

        auto fftBank = globalBankOf!(makeFFTWObject!Cpx);
        auto fftObjDn = fftBank[numOfFFT];
        auto fftObjUp = fftBank[numOfFFT * Q];
        auto dnIn = fftObjDn.inputs!F[];
        auto dnOut = fftObjDn.outputs!F[];
        auto upIn = fftObjUp.inputs!F[];
        auto upOut = fftObjUp.outputs!F[];
        auto res = new Cpx!F[][BasisFuncs.length * QX][](tx.length / (numOfFFT + numOfCp));

        foreach(r; 0 .. tx.length / (numOfFFT + numOfCp)) {
            auto sym = tx[r * (numOfFFT + numOfCp) .. $][numOfCp .. numOfCp + numOfFFT];

            dnIn[] = sym[];
            fftObjDn.fft!F();
            upIn[] = cpx!(Cpx, F)(0, 0);
            upIn[0 .. numOfFFT / 2] = dnOut[0 .. numOfFFT/2];
            upIn[$ - numOfFFT / 2 .. $] = dnOut[$ - numOfFFT / 2 .. $];
            fftObjUp.ifft!F();

            auto upSampled = fftObjUp.outputs!F.dup;

            foreach(i, BF; BasisFuncs) {
                foreach(j, e; upSampled)
                    upIn[j] = BF(e);

                // fft + swap
                fftObjUp.fft!F();
                swapHalf(fftObjUp.outputs!float);

                foreach(long j; 0 .. QX) {
                    immutable long bandIdx = (j % 2 == 0) ? j/2 : -(j+1)/2;

                    auto dst = new Cpx!F[numOfFFT];

                    if((QX == Q+1 && j < Q-1) || QX == Q)
                        dst[] = upOut[$/2 + numOfFFT * bandIdx - numOfFFT/2 .. $/2 + numOfFFT * bandIdx + numOfFFT/2];
                    else if(j == Q-1){
                        dst[] = cpx!(Cpx, F)(0, 0);
                        dst[$/2 .. $] = upOut[0 .. numOfFFT/2];
                    }else{  // if j == Q
                        dst[] = cpx!(Cpx, F)(0, 0);
                        dst[0 .. $/2] = upOut[$ - numOfFFT/2 .. $];
                    }

                    // swap + ifft
                    swapHalf(dst);
                    res[r][i*QX + j] = dst;
                }
            }
        }

        return res;
    }
}


/**
OFDM信号をオーバーサンプリングすることで折り返し雑音を再現します．
結果は時間領域です．

Q: オーバーサンプリング率
BasisFuncs: 基底関数のリスト
tx: 送信レプリカ信号, N個のシンボルを一度に変換したい場合，tx.length == N * (numOfFFT + numOfCp)
numOfFFT: OFDM信号の実効FFTサイズ
numOfCp : OFDM信号の実効CPサイズ

Return: is(typeof(result[p][i]) == typeof(tx[0]))であり，
        result[p][i]は，pバンド目の信号のi番目の時間サンプルである
*/
template generateOFDMAliasSignal(size_t Q, BasisFuncs...)
{
    import std.algorithm;
    import std.math;

    import carbon.complex;
    import dffdd.utils.fft;

    enum size_t QX = Q % 2 == 0 ? Q+1 : Q;

    C[BasisFuncs.length * QX][] generateOFDMAliasSignal(C)(in C[] tx, size_t numOfFFT, size_t numOfCp)
    if(isComplex!C)
    in{
        assert(tx.length % (numOfFFT + numOfCp) == 0);
    }
    body{
        alias Cpx = complexTypeTemplate!C;
        alias F = typeof(C.init.re);

        auto fftObj = globalBankOf!(makeFFTWObject!Cpx)[numOfFFT];

        auto freq = generateOFDMAliasSignalAtFrequency!(Q, BasisFuncs)(tx, numOfFFT, numOfCp);
        auto res = new C[typeof(freq[0]).init.length][](tx.length);

        foreach(r; 0 .. tx.length / (numOfFFT + numOfCp))
        {
            auto dst = res[r * (numOfFFT + numOfCp) .. (r+1) * (numOfFFT + numOfCp)];

            foreach(p; 0 .. freq[0].length){
                fftObj.ifftFrom!F(freq[r][p]);
                auto ops = fftObj.outputs!F;

                foreach(i; 0 .. numOfFFT)
                    dst[i + numOfCp][p] = ops[i];

                foreach(i; 0 .. numOfCp)
                    dst[i][p] = ops[$ - numOfCp + i];

              version(DFFDD_OFDM_SIDELOBE_SUPPRESS)
              {
                foreach(j, ref e; dst[$-1 .. $])
                    e[p] /= (j+2);
                foreach(j, ref e; dst[0 .. 1])
                    e[p] /= (1+1 - j);
              }
            }
        }

        return res;
    }
}


unittest
{
    import std.algorithm;
    import std.complex;
    import std.math;
    import std.range;
    import std.stdio;
    import dffdd.mod.ofdm;
    import dffdd.mod.primitives;
    import dffdd.utils.fft;


    enum int X = -100;   // そのindex, もしくはfreqに値が無いことを示す

    enum size_t numOfFFT = 8;

    /*
    |---|---|---|---|---|---|---|---|
    0   4   8   12  16  20  24  28 31 : index

    |---|---|---|---|---|---|---|---|
   16  20  24  28 31/0  4   8   12 15 : index
   -16 -12 -8  -4 -1/0  4   8   12 15 : freq

                  B   A 
      A   B   A   B   A   B   A   B 
     -2  -1  -1   0   0   1   1   2
      3     1       0        2    4
    
    */
    // 折り返しにより，indexがどこに移るか
    int[numOfFFT][5] spIdxMap4 =
        [/*         A                      B       */
            [ 0,  1,  2,  3,        28, 29, 30, 31,],
            [24, 25, 26, 27,        20, 21, 22, 23,],
            [ 8,  9, 10, 11,         4,  5,  6,  7,],
            [16, 17, 18, 19,         X,  X,  X,  X,],
            [ X,  X,  X,  X,        12, 13, 14, 15,],
        ];

    /*
    |---|---|---|---|---|---|
    0   4   8   12  16  20  24 : index

    |---|---|---|---|---|---|
   12  16  20 23/0  4   8  11 : index
   -12 -8  -4 -1/0  4   8  11 : freq

              B   A 
      B   A   B   A   B   A 
     -1  -1   0   0   1   1
        1       0        2
    
    */
    // 折り返しにより，indexがどこに移るか
    int[numOfFFT][3] spIdxMap3 =
        [/*         A                      B       */
            [ 0,  1,  2,  3,        20, 21, 22, 23,],
            [16, 17, 18, 19,        12, 13, 14, 15,],
            [ 8,  9, 10, 11,         4,  5,  6,  7,],
        ];

    static
    void tester(size_t numOfFFT, size_t numOfOS, size_t QX)(int[numOfFFT][QX] spIdxMap)
    {
        // 折り返しにより，freqがどこに移るか
        int[numOfFFT][numOfOS+1] spFrqMap;
        foreach(i; 0 .. QX) foreach(j; 0 .. numOfFFT)
            spFrqMap[i][j] = spIdxMap[i][j] >= (numOfFFT*numOfOS)/2 ? (spIdxMap[i][j] - cast(int)(numOfFFT*numOfOS)) : spIdxMap[i][j];


        auto fftObj = globalBankOf!(makeFFTWObject!(Complex))[numOfFFT];
        auto fftObjUp = globalBankOf!(makeFFTWObject!(Complex))[numOfFFT * numOfOS];

        auto inps = fftObj.inputs!float();

        foreach(A; [Complex!float(1, 0), Complex!float(SQRT1_2, SQRT1_2), Complex!float(0, 1), Complex!float(-SQRT1_2, SQRT1_2), Complex!float(-1, 0)])
            foreach(int idx; -cast(int)numOfFFT/2 .. cast(int)numOfFFT/2)
            {
                inps[] = Complex!float(0, 0);
                inps[idx < 0 ? idx + numOfFFT : idx] = A;

                fftObj.ifft!float();
                auto inpTime = fftObj.outputs!float.dup;

                // 3乗のシミュレーション結果
                fftObjUp.inputs!float[] = Complex!float(0, 0);
                fftObjUp.inputs!float[0 .. numOfFFT/2] = inps[0 .. numOfFFT/2];
                fftObjUp.inputs!float[$ - numOfFFT/2 .. $] = inps[numOfFFT/2 .. $];
                fftObjUp.ifft!float();
                fftObj.fftFrom!float(fftObjUp.outputs!float.map!(x => x*x*x).stride(numOfOS));
                auto simsig = fftObj.outputs!float.dup.map!(a => a * (numOfFFT*numOfOS)^^2 * numOfOS);

                // 3乗で試してみる
                auto aliasFreq = generateOFDMAliasSignalAtFrequency!(numOfOS, x => x*x*x)(inpTime, numOfFFT, 0);
                auto aliasTime = generateOFDMAliasSignal!(numOfOS, x => x*x*x)(inpTime, numOfFFT, 0);

                foreach(i; 0 .. QX){
                    fftObj.fftFrom!float(aliasTime.map!(a => a[i]));

                    auto ops1 = aliasFreq[0][i].map!(a => a * (numOfFFT*numOfOS)^^2);
                    auto ops2 = fftObj.outputs!float.dup.map!(a => a * (numOfFFT*numOfOS)^^2);
                    foreach(j; 0 .. numOfFFT)
                        if(spFrqMap[i][j] == idx*3){
                            auto B = (A^^3);
                            assert(approxEqual(ops1[j].re, B.re));
                            assert(approxEqual(ops1[j].im, B.im));
                            assert(approxEqual(ops2[j].re, B.re));
                            assert(approxEqual(ops2[j].im, B.im));


                            // 実際にシミュレーション結果と合致しているか確認
                            foreach(k; 0 .. numOfFFT){
                                assert(approxEqual(ops1[k].re, simsig[k].re));
                                assert(approxEqual(ops1[k].im, simsig[k].im));
                                assert(approxEqual(ops2[k].re, simsig[k].re));
                                assert(approxEqual(ops2[k].im, simsig[k].im));
                            }
                        }
                        else{
                            assert(approxEqual(ops1[j].abs(), 0));
                            assert(approxEqual(ops2[j].abs(), 0));
                        }
                }
            }
    }

    tester!(numOfFFT, 4)(spIdxMap4);
    tester!(numOfFFT, 3)(spIdxMap3);
}
