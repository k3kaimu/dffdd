module dffdd.filterbank.dftfilterbank;

import std.algorithm;
import std.array;
import std.complex;
import std.math;
import std.typecons;

import dffdd.filterbank.protofilter;
import dffdd.utils.fft;

import mir.ndslice;


/**
DFTフィルタバンクを構成します．
*/
final class DFTFilterBank(C, Flag!"isAnalysis" isAnalysis)
{
    alias R = typeof(C.init.re);


    this(size_t nchannel, in C[] prototype)
    in{
        assert(prototype.length % nchannel == 0);
    }
    do {
        immutable size_t ntaps = prototype.length / nchannel;
    
        _nchannel = nchannel;
        _ntaps = ntaps;
        _fftw = makeFFTWObject!Complex(nchannel);
        _coefs = slice!C(nchannel, ntaps);
        _inputs = slice!C(nchannel, ntaps);

        _inputs[] = C(0);
        foreach(i, e; prototype){
            immutable ci = i % nchannel;
            immutable ti = i / nchannel;

            static if(isAnalysis)
                _coefs[ci, ti] = e;
            else
                _coefs[ci, $-1-ti] = e;
        }
    }


    void opCall(in C[] inputs, ref C[] outputs)
    in {
        assert(inputs.length == _nchannel);
    }
    do {
        if(outputs.length != inputs.length)
            outputs.length = inputs.length;

        static if(isAnalysis)
        {
            // 各チャネルのFIRフィルタの状態を一つすすめる
            foreach(ic; 0 .. _nchannel)
                foreach_reverse(it; 0 .. _ntaps-1)
                    _inputs[ic, it+1] = _inputs[ic, it];

            // 入力信号を各チャネルのFIRフィルタに突っ込む
            foreach(ic; 0 .. _nchannel)
                _inputs[ic, 0] = inputs[ic];

            auto ips = _fftw.inputs!R;
            // 各チャネルで畳み込みの計算をする
            foreach(ic, ref e; _fftw.inputs!R) {
                e = 0;
                foreach(it; 0 .. _ntaps)
                    e += _inputs[ic, it] * _coefs[ic, it];
            }

            // 各チャネルの出力でIFFTをかける
            _fftw.fft!R();
            // _fftw.outputs!R[] = _fftw.inputs!R[];
            foreach(i, e; _fftw.outputs!R)
                outputs[i] = e * _nchannel;
        }
        else    // isSynthesis
        {
            // 各チャネルで入力信号をFFTする
            _fftw.inputs!R[] = inputs[];
            _fftw.ifft!R();
            // _fftw.outputs!R[] = _fftw.inputs!R[];

            // 各チャネルのFIRフィルタの状態を一つすすめる
            foreach(ic; 0 .. _nchannel)
                foreach_reverse(it; 0 .. _ntaps-1)
                    _inputs[ic, it+1] = _inputs[ic, it];

            // 入力信号を各チャネルのFIRフィルタに突っ込む
            foreach(ic, e; _fftw.outputs!R)
                _inputs[ic, 0] = e;

            // 各チャネルで畳み込みの計算をする
            foreach(ic, ref e; outputs) {
                e = 0;
                foreach(it; 0 .. _ntaps)
                    e += _inputs[ic, it] * _coefs[ic, it];
                
                e *= _nchannel;
            }
        }
    }


    size_t numOfDelay() const @property
    {
        return _nchannel * _ntaps - _nchannel;
    }


  private:
    size_t _nchannel;
    size_t _ntaps;
    FFTWObject!Complex _fftw;
    Slice!(C*, 2, Contiguous) _inputs;
    Slice!(C*, 2, Contiguous) _coefs;
}


/// ditto
alias DFTAnalysisFilterBank(C) = DFTFilterBank!(C, Yes.isAnalysis);


/// ditto
alias DFTSynthesisFilterBank(C) = DFTFilterBank!(C, No.isAnalysis);


// 
unittest
{
    alias C = Complex!double;
    auto afb = new DFTAnalysisFilterBank!C(2, [C(0.5), C(0.5)]);

    C[] inputs = new C[2],
        outputs = new C[2];
    
    inputs[0] = C(1);
    inputs[1] = C(2);

    afb(inputs, outputs);

    auto sfb = new DFTSynthesisFilterBank!C(2, [C(0.5), C(0.5)]);
    sfb(outputs, inputs);

    assert(inputs[0].re.isClose(1));
    assert(inputs[1].re.isClose(2));
    assert(inputs[0].im.isClose(0));
    assert(inputs[1].im.isClose(0));
}

// 
unittest
{
    import std.stdio;
    import std.algorithm : stdmap = map;
    import std.range;

    alias R = real;
    alias C = Complex!R;

    immutable nCH = 64,
              nTap = 32;

    auto proto = designKaiserSincPulse!R(nCH, nTap, 100).stdmap!(a => complex!R(a, 0)).array();
    // auto proto = designRootRaisedCosine!R(nCH, nTap, 0.5).stdmap!(a => complex!R(a, 0)).array();
    // auto proto = designRootRaisedERF!R(nCH/2, nTap*2, R.nan, 0.911).stdmap!(a => complex!R(a, 0)).array();
    assert(proto.length == nCH * nTap);
    immutable delay = nCH * nTap - nCH;
    immutable len = nCH*nTap*2;

    auto afb = new DFTAnalysisFilterBank!C(nCH, proto);
    auto sfb = new DFTSynthesisFilterBank!C(nCH, proto);

    auto chirp = sequence!"(a[0]*n*(n-1)/2) % 1"(1.0/len).stdmap!(a => std.complex.expi(-a * 2 * PI)).take(len).array();
    C[] b1 = new C[nCH],
        b2 = new C[nCH];
    
    C[] dst = new C[len];
    foreach(i; 0 .. nTap*2) {
        b1[0 .. nCH] = chirp[i*nCH .. (i+1)*nCH];
        afb(b1, b2);
        sfb(b2, b1);
        dst[i*nCH .. (i+1)*nCH] = b1[];
    }

    dst = dst[delay .. $];
    chirp = chirp[0 .. dst.length];

    real sum = 0, P = 0;
    foreach(e; dst.zip(chirp)){
        sum += (e[0] - e[1]).sqAbs;
        P += e[0].sqAbs;
    }

    // reconstruction error is less than -40 dB.
    assert(10*log10(sum / P) < -40);
}


/**
See also: https://jp.mathworks.com/matlabcentral/fileexchange/15813-near-perfect-reconstruction-polyphase-filterbank
*/
final class OversampledDFTFilterBank(C, Flag!"isAnalysis" isAnalysis)
{
    this(size_t nchannel, in C[] prototype)
    {
        _nchannel = nchannel;
        _bank0 = new DFTFilterBank!(C, isAnalysis)(nchannel/2, prototype);
        _bank1 = new DFTFilterBank!(C, isAnalysis)(nchannel/2, prototype);

        _oddCoefs = new C[_nchannel/2];
        foreach(i, ref e; _oddCoefs)
            e = std.complex.expi(PI * i / (_nchannel/2));

        _buf0 = new C[_nchannel/2];
        _buf1 = new C[_nchannel/2];
        _buf2 = new C[_nchannel/2];
    }


    void opCall(in C[] inputs, ref C[] outputs)
    in {
        if(isAnalysis){
            assert(inputs.length == _nchannel / 2);
        }else{
            assert(inputs.length == _nchannel);
        }
    }
    do {
        if(isAnalysis) {
            outputs.length = _nchannel;
        } else {
            outputs.length = _nchannel / 2;
        }

        static if(isAnalysis)
        {
            // even channel
            _bank0(inputs, _buf0);
            foreach(i, e; _buf0)
                outputs[i*2] = e;

            // odd channel
            foreach(i, ref e; _buf1){
                e = inputs[i] * _oddCoefs[i];
                if(_opCallCNT % 2 == 1)
                    e *= -1;
            }

            ++_opCallCNT;
            _bank1(_buf1, _buf0);
            foreach(ptrdiff_t i, e; _buf0)
                outputs[((i*2-1) + $) % $] = e;
        }
        else
        {
            // split to odd and even elements
            foreach(ptrdiff_t i; 0 .. _nchannel / 2){
                _buf0[i] = inputs[i*2];
                _buf1[i] = inputs[((i*2-1) + $) % $ ];
            }

            // even channel
            _bank0(_buf0, outputs);

            // odd channel
            _bank1(_buf1, _buf0);
            foreach(i, ref e; _buf0){
                e = _buf0[i] * _oddCoefs[i].conj();

                if(_opCallCNT % 2 == 1)
                    e *= -1;
            }

            ++_opCallCNT;

            // merge two filter banks
            foreach(i, ref e; outputs)
                e -= _buf0[i];
        }
    }


    size_t numOfDelay() const @property
    {
        return _bank0.numOfDelay;
    }


    size_t numOfChannel() const @property
    {
        return this._nchannel;
    }


  private:
    DFTFilterBank!(C, isAnalysis) _bank0, _bank1;
    size_t _nchannel;
    size_t _opCallCNT;
    C[] _oddCoefs;
    C[] _buf0, _buf1, _buf2;
}


/// ditto
alias OversampledDFTAnalysisFilterBank(C) = OversampledDFTFilterBank!(C, Yes.isAnalysis);
alias OversampledDFTSynthesisFilterBank(C) = OversampledDFTFilterBank!(C, No.isAnalysis);


// 
unittest
{
    import std.stdio;
    import std.algorithm : stdmap = map;
    import std.range;
    import dffdd.window;

    alias R = real;
    alias C = Complex!R;

    immutable nCH = 64,
              nTap = 32;

    auto proto = designKaiserBesselDrived!R(nCH, nTap/2, kaiserBetaFromStopbandAttdB(200), 1).stdmap!(a => complex!R(a, 0)).array();
    // auto proto = designRootRaisedERF!R(nCH/2, nTap).stdmap!(a => complex!R(a, 0)).array();
    // auto proto = designKaiserSincPulse!R(nCH/2, nTap, 50, 1.9125).stdmap!(a => complex!R(a, 0)).array();
    // auto proto = designRootRaisedCosine!R(nCH, nTap/2, 0.5).stdmap!(a => complex!R(a, 0)).array();
    assert(proto.length == nCH/2 * nTap);
    immutable delay = nCH/2 * nTap - nCH/2;
    immutable len = nCH*nTap;

    auto afb = new OversampledDFTAnalysisFilterBank!C(nCH, proto);
    auto sfb = new OversampledDFTSynthesisFilterBank!C(nCH, proto);

    auto chirp = sequence!"(a[0]*n*(n-1)/2) % 1"(1.0/len).stdmap!(a => std.complex.expi(-a * 2 * PI)).take(len).array();
    C[] b1 = new C[nCH/2],
        b2 = new C[nCH];

    C[] dst = new C[len];
    foreach(i; 0 .. nTap*2) {
        b1[0 .. nCH/2] = chirp[i*nCH/2 .. (i+1)*nCH/2];
        afb(b1, b2);
        sfb(b2, b1);
        dst[i*nCH/2 .. (i+1)*nCH/2] = b1[];
    }

    dst = dst[delay .. $];
    chirp = chirp[0 .. dst.length];

    real sum = 0, P = 0;
    foreach(e; dst.zip(chirp)){
        sum += (e[0] - e[1]).sqAbs;
        P += e[0].sqAbs;
    }

    // writeln(10*log10(sum / P));
    assert(10*log10(sum / P) < -150);
}


final class OversampledDFTFilterBankModulator(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(size_t nchannel, in C[] prototype)
    {
        this._analysis = new OversampledDFTAnalysisFilterBank!C(nchannel, prototype);
        this._synthesis = new OversampledDFTSynthesisFilterBank!C(nchannel, prototype);
    }


    size_t symInputLength() const @property
    {
        return this._analysis.numOfChannel();
    }


    size_t symOutputLength() const @property
    {
        return this._analysis.numOfChannel() / 2;
    }


    ref OutputElementType[] modulate(in InputElementType[] inputs, return ref OutputElementType[] outputs)
    in {
        assert(inputs.length % this.symInputLength == 0);
    }
    do {
        if(outputs.length != inputs.length / this.symInputLength * this.symOutputLength)
            outputs.length = inputs.length / 2;
        
        foreach(i; 0 .. _inputs.length / this.symInputLength) {
            _analysis(
                inputs[i * this.symInputLength .. (i+1) * this.symInputLength],
                outputs[i * this.symOutputLength .. (i+1) * this.symOutputLength]);
        }
    }

    ref InputElementType[] demodulate(in OutputElementType[] inputs, return ref InputElementType[] outputs)
    in {
        assert(inputs.length % this.symOutputLength == 0);
    }
    do {
        if(outputs.length != inputs.length / this.symOutputLength * this.symInputLength)
            outputs.length = inputs.length / this.symOutputLength * this.symInputLength;

        foreach(i; 0 .. _inputs.length / this.symOutputLength) {
            _synthesis(
                inputs[i * this.symOutputLength .. (i+1) * this.symOutputLength],
                outputs[i * this.symInputLength .. (i+1) * this.symInputLength]
            );
        }
    }


  private:
    OversampledDFTAnalysisFilterBank!C _analysis;
    OversampledDFTSynthesisFilterBank!C _synthesis;
}

// unittest
// {
//     import std.algorithm;
//     import std.format;
//     import std.stdio;

//     auto fftw = makeFFTWObject!Complex(1024*16);

//     void outputIRAndFR(real[] coefs, string filename)
//     {
//         File(filename ~ "_time.csv", "w").writefln!"[%(%.20g,%)]"(coefs);
//         fftw.inputs!real[] = Complex!real(0);
//         foreach(i, ref e; fftw.inputs!real[0 .. coefs.length]) e = coefs[i];
//         fftw.fft!real();
//         File(filename ~ "_freq.csv", "w").writefln!"[%(%.20g,%)]"(fftw.outputs!real.map!sqAbs);
//     }

//     foreach(beta; [0, 0.25, 0.5, 0.75, 1]) {
//         auto coefs = designRootRaisedCosine!real(32, 32, beta);
//         outputIRAndFR(coefs, "rrc_%s".format(beta));
//     }

//     foreach(beta; [2, 4, 8, 16, 32]) {
//         auto coefs = designKaiserBesselDrived!real(32, 32, beta);
//         outputIRAndFR(coefs, "kbd_%s".format(beta));
//     }

//     foreach(k; [4, 6, 8, 10, 12]) {
//         auto coefs = designRootRaisedERF!real(32/2, 64, k);
//         outputIRAndFR(coefs, "erf_%s".format(k));
//     }
// }