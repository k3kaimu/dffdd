module dffdd.mod.sefdm;

import std.exception : assumeUnique;
import std.random;
import std.traits;


import dffdd.mod.primitives : Bit;
import dffdd.mod.ofdm;
import dffdd.math.complex;
import dffdd.math.math;
import dffdd.utils.fft;


final class RDFTsSEFDM(Mod, C = Mod.OutputELementType)
{
    alias InputElementType = Bit;
    alias OutputElementType = C;
    alias RealType = typeof(C.init.re);


    this(
        Mod mod, uint nFFT, uint nCP, uint nTone, uint nUpSampling,
        uint nData, immutable(size_t)[] shuffle = null, uint shuffleRandSeed = 0,
        double N0 = 1e-6, size_t maxIter = 20)
    in(shuffle.length == 0 || shuffle.length == nData)
    in(nData >= nTone)
    in(nFFT >= nTone)
    {
        _mod = mod;
        _ofdm = new OFDM!C(nFFT, nCP, nTone, nUpSampling);
        _nTone = nTone;
        _nData = nData;
        _cbuffer = new C[_nData];
        _sbuffer = new C[_nTone];

        _fftw = makeFFTWObject!C(_nData);

        // シャッフルするインデックス配列の設定
        if(shuffle !is null) {
            _rndIdx = shuffle;
        } else {
            // 引数で与えられていない場合は自分で生成する
            auto perm = new size_t[nData];
            foreach(i; 0 .. nData)
                perm[i] = i;
            
            Random rnd;
            rnd.seed(shuffleRandSeed);
            perm = randomShuffle(perm, rnd);
            _rndIdx = assumeUnique(perm);
        }

        _epdet = _makeEPDet(N0, maxIter);
    }


    size_t symInputLength() const @property { return _nData * _mod.symInputLength; }
    size_t symOutputLength() const @property { return _ofdm.symOutputLength; }


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


  private:
    alias F = typeof(C.init.re);

    Mod _mod;
    OFDM!C _ofdm;
    typeof(makeFFTWObject!C(0)) _fftw;
    C[] _sbuffer, _cbuffer;
    size_t _nTone;
    size_t _nData;
    immutable(size_t)[] _rndIdx;
    typeof(_makeEPDet(0, 0)) _epdet;


    auto _makeEPDet(double N0, size_t maxIter)
    {
        import dffdd.math.matrix;
        import dffdd.math.vector;
        import dffdd.math.matrixspecial : PermutationMatrix, dftMatrix, identity;
        import dffdd.detector.ep : makeSVDEPDetector;

        auto pmat = PermutationMatrix(_rndIdx);
        auto wmat = dftMatrix!C(_nData);
        return makeSVDEPDetector!C(_mod, identity!C(_nTone), vector!C(_nTone, C(1)), pmat * wmat, N0, maxIter);
    }


    void _randomSpread(in C[] inbuf, C[] outbuf)
    in(inbuf.length == _nData)
    in(outbuf.length == _nTone)
    {
        alias F = typeof(C.init.re);
        immutable F scale = 1/fast_sqrt!F(inbuf.length);

        _fftw.inputs!F[] = inbuf[];
        _fftw.fft!F();
        foreach(i; 0 .. _nTone)
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


    void runTestNoDist(uint nFFT, uint nCP, uint nTone, uint nUpSampling, uint nData, uint nSym)
    {
        auto sefdm = new RDFTsSEFDM!(QPSK!C, C)(qpsk, nFFT, nCP, nTone, nUpSampling, nData, null, 0);
        foreach(_; 0 .. 100) {
            auto bits = randomBits().map!(a => Bit(a)).take(nData * 2 * nSym).array();

            C[] dst_mod;
            sefdm.modulate(bits, dst_mod);
            assert(dst_mod.length == (nFFT+nCP)*nSym);

            Bit[] dst_demod;
            sefdm.demodulate(dst_mod, dst_demod);
            assert(dst_demod.length == nData * 2 * nSym);
            assert(bits == dst_demod);
        }
    }


    runTestNoDist(64, 16, 52, 1, 80, 1);
    runTestNoDist(64, 0, 64, 1, 64, 2);
    runTestNoDist(2, 0, 2, 1, 2, 100);


    C[] toLTI(in C[] signal, in C[] taps)
    {
        C[] recv = new C[signal.length];
        recv[] = C(0);
        foreach(n, e; taps)
            foreach(i; n .. recv.length)
                recv[i] += e * signal[i - n];
        
        return recv;
    }


    void runTestLTI(uint nFFT, uint nCP, uint nTone, uint nUpSampling, uint nData, uint nSym, in C[] taps)
    {
        auto sefdm = new RDFTsSEFDM!(QPSK!C, C)(qpsk, nFFT, nCP, nTone, nUpSampling, nData, null, 0);

        // チャネルの周波数応答を推定する
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

            sefdm.epDetector.channelSingularValues[] = ch.sliced.vectored;
        }

        foreach(_; 0 .. 100) {
            auto bits = randomBits().map!(a => Bit(a)).take(nData * 2 * nSym).array();

            C[] dst_mod;
            sefdm.modulate(bits, dst_mod);
            assert(dst_mod.length == (nFFT+nCP)*nSym);

            C[] recv = toLTI(dst_mod, taps);

            Bit[] dst_demod;
            sefdm.demodulate(recv, dst_demod);
            assert(dst_demod.length == nData * 2 * nSym);
            assert(bits == dst_demod);
        }
    }


    runTestLTI(64, 16, 52, 1, 52, 1, [C(1), C(0.1), C(-0.1)]);
    runTestLTI(64, 16, 52, 1, 80, 1, [C(0.1), C(1), C(-0.1), C(0.1), C(-0.1), C(0.1), C(-0.1)]);
}