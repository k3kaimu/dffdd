module dffdd.mod.sefdm;

import std.exception : assumeUnique;
import std.random;
import std.traits;


import dffdd.mod.primitives : Bit;
import dffdd.mod.ofdm;
import dffdd.math.complex;
import dffdd.math.math;
import dffdd.utils.fft;


final class RDFTsSEFDM(Mod, C = Mod.OutputElementType)
{
    alias InputElementType = Bit;
    alias OutputElementType = C;
    alias RealType = typeof(C.init.re);


    this(
        Mod mod, uint nFFT, uint nCP, uint nTone, uint nUpSampling,
        uint nSpreadIn, uint nSpreadOut, immutable(size_t)[] shuffle = null, uint shuffleRandSeed = 0,
        double N0 = 1e-6, size_t maxIter = 20)
    in(shuffle.length == 0 || shuffle.length == nSpreadIn)
    in(nUpSampling >= 1)
    in(nFFT >= nTone)
    in(nSpreadIn >= nSpreadOut)
    in(nSpreadOut % nTone == 0)
    {
        _mod = mod;
        _ofdm = new OFDM!C(nFFT, nCP, nTone, nUpSampling);
        _nTone = nTone;
        _nFFT = nFFT;
        _nSpreadIn = nSpreadIn;
        _nSpreadOut = nSpreadOut;
        _cbuffer = new C[_nSpreadIn];
        _sbuffer = new C[_nSpreadOut];

        _fftw = makeFFTWObject!C(_nSpreadIn);

        // シャッフルするインデックス配列の設定
        if(shuffle !is null) {
            _rndIdx = shuffle;
        } else {
            // 引数で与えられていない場合は自分で生成する
            auto perm = new size_t[nSpreadIn];
            foreach(i; 0 .. nSpreadIn)
                perm[i] = i;
            
            Random rnd;
            rnd.seed(shuffleRandSeed);
            perm = randomShuffle(perm, rnd);
            _rndIdx = assumeUnique(perm);
        }

        // OFDM.demodulateの中で係数がかけられるのでその分だけN0に補正をかけておく
        _epdet = _makeEPDet(N0 * nTone / nFFT, maxIter);
    }


    size_t symInputLength() const @property
    {
        return _mod.symInputLength * _nSpreadIn;
    }


    size_t symOutputLength() const @property
    {
        immutable nblock = _nSpreadOut / _nTone;
        return _ofdm.symOutputLength * nblock;
    }


    ref C[] modulate(in Bit[] inputs, return ref C[] outputs)
    in(inputs.length % this.symInputLength == 0)
    do {
        immutable M = this.symInputLength;
        immutable N = this.symOutputLength;

        if(outputs.length != inputs.length / M * N)
            outputs.length = inputs.length / M * N;

        foreach(i; 0 .. inputs.length / M) {
            auto inbuf = inputs[i * M .. (i+1) * M];
            auto outbuf = outputs[i * N .. (i+1) * N];

            _cbuffer = _mod.modulate(inbuf, _cbuffer);
            _randomSpread(_cbuffer, _sbuffer);
            _ofdm.modulate(_sbuffer, outbuf);
        }

        return outputs;
    }


    ref C[] demodulateAsOFDM(in C[] received, return ref C[] subcarriers)
    in(received.length % this.symOutputLength == 0)
    {
        return _ofdm.demodulate(received, subcarriers);      
    }


    ref Bit[] demodulate(in C[] received, return ref Bit[] detected)
    in(received.length % this.symOutputLength == 0)
    {
        immutable M = this.symInputLength;
        immutable N = this.symOutputLength;

        if(detected.length != received.length / N * M)
            detected.length = received.length / N * M;

        foreach(n; 0 .. received.length / N) {
            // OFDMとして復調
            _ofdm.demodulate(received[n*N .. (n+1)*N], _sbuffer);
            // EPで復調
            _epdet.detect(_sbuffer, _cbuffer);
            // サブキャリアの復調
            auto dst1 = detected[n*M .. (n+1)*M];
            _mod.demodulate(_cbuffer, dst1);
        }

        return detected;
    }


    void setFrequencyResponse(in C[] freqResp)
    in(freqResp.length == _nTone)
    {
        import mir.ndslice;
        auto fr = freqResp.sliced;

        foreach(i; 0 .. _nSpreadOut / _nTone)
            _epdet.channelSingularValues.sliced()[i * _nTone .. (i+1) * _nTone] = freqResp.sliced;
    }


    void setNoisePower(double N0)
    {
        _epdet.setNoisePower(N0 * _nTone / _nFFT);
    }


    ref inout(OFDM!C) ofdm() inout
    {
        return _ofdm;
    }


    ref inout(Mod) scmod() inout
    {
        return _mod;
    }


    ref inout(typeof(_epdet)) epDetector() inout
    {
        return _epdet;
    }


    immutable(size_t)[] scShuffle() const { return _rndIdx; }


  private:
    alias F = typeof(C.init.re);

    Mod _mod;
    OFDM!C _ofdm;
    typeof(makeFFTWObject!C(0)) _fftw;
    C[] _sbuffer, _cbuffer;
    size_t _nTone, _nFFT;
    size_t _nSpreadIn, _nSpreadOut;
    immutable(size_t)[] _rndIdx;
    typeof(_makeEPDet(0, 0)) _epdet;


    auto _makeEPDet(double N0, size_t maxIter)
    {
        import dffdd.math.matrix;
        import dffdd.math.vector;
        import dffdd.math.matrixspecial : PermutationMatrix, dftMatrix, identity;
        import dffdd.detector.ep : makeSVDEPDetector;

        auto pmat = PermutationMatrix(_rndIdx);
        auto wmat = dftMatrix!C(_nSpreadIn);
        return makeSVDEPDetector!C(_mod, identity!C(_nSpreadOut), vector!C(_nSpreadOut, C(1)), pmat * wmat, N0, maxIter);
    }


    void _randomSpread(in C[] inbuf, C[] outbuf)
    in(inbuf.length == _nSpreadIn)
    in(outbuf.length == _nSpreadOut)
    {
        alias F = typeof(C.init.re);
        immutable F scale = 1/fast_sqrt!F(_nSpreadIn);

        _fftw.inputs!F[] = inbuf[];
        _fftw.fft!F();
        foreach(i; 0 .. _nSpreadOut)
            outbuf[i] = _fftw.outputs!F[_rndIdx[i]] * scale;
    }
}


unittest
{
    import std.algorithm;
    import std.range;
    import dffdd.utils.binary;
    import dffdd.mod.qpsk;

    alias C = Complex!float;
    alias T = RDFTsSEFDM!(QPSK!C, C);
    QPSK!C qpsk;


    void runTestNoDist(uint nFFT, uint nCP, uint nTone, uint nUpSampling, uint nSpreadIn, uint nSpreadOut, uint nSym)
    {
        auto sefdm = new RDFTsSEFDM!(QPSK!C, C)(qpsk, nFFT, nCP, nTone, nUpSampling, nSpreadIn, nSpreadOut, null, 0);
        foreach(_; 0 .. 100) {
            auto bits = randomBits().map!(a => Bit(a)).take(nSpreadIn * 2 * nSym).array();

            C[] dst_mod;
            sefdm.modulate(bits, dst_mod);
            assert(dst_mod.length == (nFFT+nCP) * nUpSampling * (nSpreadOut / nTone) * nSym);

            Bit[] dst_demod;
            sefdm.demodulate(dst_mod, dst_demod);
            assert(dst_demod.length == nSpreadIn * 2 * nSym);
            assert(bits == dst_demod);
        }
    }


    runTestNoDist(64, 16, 52, 4, 53, 52, 100);
    runTestNoDist(64, 0, 64, 1, 64, 64, 100);
    runTestNoDist(2, 0, 2, 1, 2, 2, 100);
    runTestNoDist(64, 16, 52, 4, 800, 520, 1);


    C[] toLTI(in C[] signal, in C[] taps)
    {
        C[] recv = new C[signal.length];
        recv[] = C(0);
        foreach(n, e; taps)
            foreach(i; n .. recv.length)
                recv[i] += e * signal[i - n];
        
        return recv;
    }


    void runTestLTI(uint nFFT, uint nCP, uint nTone, uint nUpSampling, uint nSpreadIn, uint nSpreadOut, uint nSym, in C[] taps)
    {
        auto sefdm = new RDFTsSEFDM!(QPSK!C, C)(qpsk, nFFT, nCP, nTone, nUpSampling, nSpreadIn, nSpreadOut, null, 0);

        // 1000個のOFDMシンボルでチャネルの周波数応答を推定する
        {
            import dffdd.mod.primitives : mod, demod;
            import mir.ndslice;
            import dffdd.math.vector;
            auto bits = randomBits().map!(a => Bit(a)).take(nTone * 2 * 1000).array();
            auto scstx = mod(sefdm.scmod, bits);
            auto ofdmtx = mod(sefdm.ofdm, scstx);
            auto ofdmrx = toLTI(ofdmtx, taps);
            auto scsrx = demod(sefdm.ofdm, ofdmrx);

            foreach(i; 0 .. scsrx.length)
                scsrx[i] /= scstx[i];

            C[] ch = new C[nTone];
            ch[] = C(0);
            foreach(i; 0 .. nTone)
                foreach(j; 0 .. scsrx.length / nTone)
                    ch[i] += scsrx[j*nTone + i] / (scsrx.length / nTone);

            sefdm.setFrequencyResponse(ch);
        }

        foreach(_; 0 .. 100) {
            auto bits = randomBits().map!(a => Bit(a)).take(nSpreadIn * 2 * nSym).array();

            C[] dst_mod;
            sefdm.modulate(bits, dst_mod);
            assert(dst_mod.length == (nFFT+nCP) * nUpSampling * (nSpreadOut / nTone) * nSym);

            C[] recv = toLTI(dst_mod, taps);

            Bit[] dst_demod;
            sefdm.demodulate(recv, dst_demod);
            assert(dst_demod.length == nSpreadIn * 2 * nSym);
            assert(bits == dst_demod);
        }
    }


    runTestLTI(64, 16, 52, 1, 52, 52, 100, [C(1), C(0.1), C(-0.1)]);
    runTestLTI(64, 16, 52, 4, 53, 52, 100, [C(0.1), C(1), C(-0.1), C(0.1), C(-0.1), C(0.1), C(-0.1)]);
    runTestLTI(64, 16, 52, 1, 800, 520, 1, [C(0.1), C(1), C(-0.1), C(0.1), C(-0.1), C(0.1), C(-0.1)]);


    // auto sefdm = new RDFTsSEFDM!(QPSK!C, C)(qpsk, 32, 0, 20, 1, 20, null, 0);
    // C[] xs = [C(1)] ~ iota(31).map!(a => C(0)).array();
    // C[] ys;
    // sefdm.ofdm.demodulate(xs, ys);
    // import std.stdio;
    // writeln(ys);
}