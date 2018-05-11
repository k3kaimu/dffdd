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


    void opCall(in C[] inputs, C[] outputs)
    in {
        assert(inputs.length == outputs.length);
        assert(inputs.length == _nchannel);
    }
    do {
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
            _fftw.ifft!R();
            foreach(i, e; _fftw.outputs!R)
                outputs[i] = e * _nchannel;
        }
        else    // isSynthesis
        {
            // 各チャネルで入力信号をFFTする
            _fftw.inputs!R[] = inputs[];
            _fftw.fft!R();

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


  private:
    size_t _nchannel;
    size_t _ntaps;
    FFTWObject!Complex _fftw;
    ContiguousSlice!(2, C) _inputs;
    ContiguousSlice!(2, C) _coefs;
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

    assert(inputs[0].re.approxEqual(1));
    assert(inputs[1].re.approxEqual(2));
    assert(inputs[0].im.approxEqual(0));
    assert(inputs[1].im.approxEqual(0));
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
    // writeln(10*log10(sum / P));
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


    void opCall(in C[] inputs, C[] outputs)
    in {
        if(isAnalysis){
            assert(inputs.length == _nchannel / 2);
            assert(outputs.length == _nchannel);
        }else{
            assert(inputs.length == _nchannel);
            assert(outputs.length == _nchannel / 2);
        }
    }
    do {
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
            foreach(i, e; _buf0)
                outputs[i*2+1] = e;
        }
        else
        {
            // split to odd and even elements
            foreach(i; 0 .. _nchannel / 2){
                _buf0[i] = inputs[i*2];
                _buf1[i] = inputs[i*2+1];
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

    alias R = real;
    alias C = Complex!R;

    immutable nCH = 64,
              nTap = 32;

    auto proto = designRootRaisedERF!R(nCH/2, nTap).stdmap!(a => complex!R(a, 0)).array();
    // auto proto = designKaiserSinc!R(nCH/2, nTap, 50, 1.9125).stdmap!(a => complex!R(a, 0)).array();
    // auto proto = designRootRaisedCosine!R(nCH, nTap/2, 1).stdmap!(a => complex!R(a, 0)).array();
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

    assert(10*log10(sum / P) < -130);
}





__EOF__


auto mapSlice(alias fn, X)(X x)
{
    return mir.ndslice.topology.map!fn(x);
}


final class NPROversampledDFTFilterBank(C, alias makeFFTObj)
{
    alias InputElementType = C;
    alias OutputElementType = C;

    alias F = typeof(C.init.re);

    this(size_t nchannel, size_t noverlap)
    {
        _setParam(nchannel, noverlap);
    }


    size_t nchannel() const @property { return _nchannel; }
    size_t noverlap() const @property { return _noverlap; }
    size_t delay() const @property { return (_nchannel / 2) * (_noverlap - 1) / 2; }

    size_t symInputLength() const @property { return this._nchannel / 2; }
    size_t symOutputLength() const @property { return this._nchannel; }


    /**
    C[N] -> C[2N]
    */
    C[] analysis(in C[] input, return ref C[] output)
    in{
        // static if(hasLength!R)
        //     assert(input.length % this.nchannel == 0);
        assert(input.length % (this.nchannel/2) == 0);
        assert(output.length % this.nchannel == 0);
    }
    body{
        immutable N = this.nchannel / 2,
                  M = input.length / N;

        Slice!(2, Unqual!C*) x1 = input.sliced(N*M).matlabReshape(N, M).slice;
        Slice!(2, Unqual!C*) x2 = slice!(Unqual!C)(N, M);
        x2[] = x1;

        foreach(i; 0 .. N){
            x2[i][] = x2[i].mapSlice!(a => a * std.complex.expi(PI * i / N));
            x2[i, 1 .. $].strided!0(2)[] = x2[i, 1 .. $].strided!0(2).mapSlice!"-a";
        }

        auto c1 = coeff.reversed!1;
        foreach(i; 0 .. N){
            x1[i][] = x1[i].firFiltered(c1[i]);
            x2[i][] = x2[i].firFiltered(c1[i]);
        }

        foreach(j; 0 .. M){
            .ifft!F(globalBankOf!makeFFTObj[N], x1[0 .. $, j], x1[0 .. $, j]);
            .ifft!F(globalBankOf!makeFFTObj[N], x2[0 .. $, j], x2[0 .. $, j]);
            x1[0 .. $, j] *= N;
            x2[0 .. $, j] *= N;
        }

        Slice!(2, Unqual!C*) ret = slice!(Unqual!C)(N*2, M);

        ret.strided!0(2)[] = x1;
        ret[1 .. $].strided!0(2)[] = x2;

        return ret;
    }


    alias demodulate = analysis;
    alias modulate = synthesis;


  private:
    Slice!(2, C*) _coeff;
    size_t _nchannel, _noverlap;
    C[][] _firStateA, _firStateS;


    void setParam(size_t N, size_t L)
    {
        _nchannel = N * 2;
        _noverlap = L;
        _coeff = calculateCoeff!F(N, L);
    }


    void analysisImpl(in C[] input, C[] output)
    in{
        assert(input.length == this.symInputLength);
        assert(output.length == this.symOutputLength);
    }
    body{
        
    }
}


auto firFiltered(C, F, SliceKind kind1, SliceKind kind2)(Slice!(kind1, [1], C*) input, Slice!(kind2, [1], F*) coeff)
{
    auto state = new Unqual!C[](coeff.length);
    auto output = slice!(Unqual!C)(input.length);
    state[0 .. $] = C(0, 0);
    foreach(i; 0 .. input.length){
        foreach_reverse(j; 0 .. state.length - 1)
            state[j+1] = state[j];
        
        state[0] = input[i];
        output[i] = 0;
        foreach(j; 0 .. state.length)
            output[i] += state[j] * coeff[j];
    }

    return output;
}


auto analysis(alias makeFFTObj = makeFFTWObject, C, F, SliceKind kind)(Slice!(kind, [2], F*) coeff, in C[] input)
in{
    assert(input.length % coeff.length!0 == 0);
}
body{
    immutable N = coeff.length!0,
              M = input.length / N;

    auto x1 = input.dup.sliced(N*M).matlabReshape(N, M).slice;
    auto x2 = slice!(Unqual!C)(N, M);
    x2[] = x1;

    foreach(i; 0 .. N){
        immutable c = std.complex.expi(PI * i / N);
        foreach(j; 0 .. x2[i].length)
            x2[i, j] *= c;
        
        foreach(j; 0 .. x2[i].length)
            if(j % 2 == 1)
                x2[i, j] *= -1;
    }

    //writeln(x2);

    auto c1 = coeff.universal.reversed!1;
    foreach(i; 0 .. N){
        x1[i, 0 .. $] = x1[i].firFiltered(c1[i]);
        x2[i, 0 .. $] = x2[i].firFiltered(c1[i]);
    }

    foreach(j; 0 .. M){
        inplaceIFFT!F(globalBankOf!makeFFTObj[N], x1[0 .. $, j]);
        inplaceIFFT!F(globalBankOf!makeFFTObj[N], x2[0 .. $, j]);
        x1[0 .. $, j] *= N;
        x2[0 .. $, j] *= N;
    }

    auto ret = slice!(Unqual!C)(N*2, M).universal;

    ret.strided!0(2)[] = x1;
    ret[1 .. $].strided!0(2)[] = x2;

    return ret;
}


void inplaceIFFT(F, FFTObject, Slice)(FFTObject fft, Slice slice)
{
    foreach(i, ref e; fft.inputs!F)
        e = slice[i];

    fft.ifft!F();

    foreach(i, e; fft.outputs!F)
        slice[i] = e;
}


auto synthesis(alias makeFFTObj = makeFFTWObject, C, R, SliceKind kind1, SliceKind kind2)(Slice!(kind1, [2], R*) coeff, Slice!(kind2, [2], C*) input)
{
    immutable N = coeff.length!0,
              M = input.length!1;

    auto y1 = input.strided!0(2).slice.universal,
         y2 = input[1 .. $].strided!0(2).slice.universal;

    foreach(j; 0 .. M){
        .fft!R(globalBankOf!makeFFTObj[N], y1[0 .. $, j], y1[0 .. $, j]);
        .fft!R(globalBankOf!makeFFTObj[N], y2[0 .. $, j], y2[0 .. $, j]);
        y1[0 .. $, j] *= N;
        y2[0 .. $, j] *= N;
    }
    //writeln(y1[0 .. 4, 0 .. 4]);

    foreach(i; 0 .. N){
        y1[i][] = y1[i].firFiltered(coeff[i]);
        y2[i][] = y2[i].firFiltered(coeff[i]);
    }

    foreach(i; 0 .. N){
        y2[i][] = y2[i].mapSlice!(a => a * std.complex.expi(- PI * i / N));
        y2[i, 1 .. $].strided!0(2)[] = y2[i, 1 .. $].strided!0(2).mapSlice!"-a";
    }

    y1[] -= y2;
    return y1.transposed.flattened;
}


private
auto rrerf(S, F)(S freq, F K, size_t M)
{
    alias mmap = mir.ndslice.topology.map;

    auto ns = new Complex!(ElementType!(S))[freq.length];
    foreach(i; 0 .. freq.length)
        ns[i] = complex!F(sqrt(erfc(K*(2*M*freq[i] - 0.5))*0.5), 0);

    return ns.sliced;
}


auto calculateCoeff(F = real, alias makeFFTObj = makeFFTWObject)(size_t N = 256, size_t L = 128, F K = F.nan)
{
    if(K.isNaN){
        switch(L){
            case 8:     K = 4.853; break;
            case 10:    K = 4.775; break;
            case 12:    K = 5.257; break;
            case 14:    K = 5.736; break;
            case 16:    K = 5.856; break;
            case 18:    K = 7.037; break;
            case 20:    K = 6.499; break;
            case 22:    K = 6.483; break;
            case 24:    K = 7.410; break;
            case 26:    K = 7.022; break;
            case 28:    K = 7.097; break;
            case 30:    K = 7.755; break;
            case 32:    K = 7.452; break;
            case 48:    K = 8.522; break;
            case 64:    K = 9.396; break;
            case 96:    K = 10.785; break;
            case 128:   K = 11.5; break; 
            case 192:   K = 11.5; break;
            case 256:   K = 11.5; break;
            default:    K = 8;
        }
    }

    immutable M = N / 2;
    auto freq = mir.ndslice.topology.iota(L*M).unaryFun!(a => mir.ndslice.topology.map!(a => cast(F)a / (L*M))(a));
    auto a = rrerf(freq, K, M).universal;
    a[$/2 .. $].reversed!0[] = mir.ndslice.topology.map!(a => a.conj)(a[1 .. $/2+1]);
    a[$/2] = 0;
    //writeln("START----------");
    //writefln("%(%s,%|\n%)", a.map!"a.re^^2 + a.im^^2");
    //writeln("END------------");
    .ifft!F(globalBankOf!makeFFTObj[a.length], a, a);
    swapHalf(a);
    a[] /= sum(a);

    foreach(e; a)
        assert(approxEqual(e.im, 0));

    return mir.ndslice.topology.map!"a.re"(a.matlabReshape(M, L)).slice;
}

unittest
{
    import std.stdio;

    auto fftobj64 = makeFFTWObject!Complex(64);
    auto fftobj32 = makeFFTWObject!Complex(32);
    auto coeff = calculateCoeff(64, 32);

    auto testResult = 
        [2.26E-09, 4.42E-09, 5.90E-09, -1.20E-08, -7.56E-08, -2.05E-07, -2.37E-07, 2.74E-06, 1.11E-05, -2.49E-05, -0.000111495, 0.000219418, 0.000564861, -0.001469149, -0.001504252, 0.009080523, 0.017713452, 0.009080523, -0.001504252, -0.001469149, 0.000564861, 0.000219418, -0.000111495, -2.49E-05, 1.11E-05, 2.74E-06, -2.37E-07, -2.05E-07, -7.56E-08, -1.20E-08, 5.90E-09, 4.42E-09,
        2.27E-09, 4.53E-09, 5.77E-09, -1.32E-08, -7.85E-08, -2.11E-07, -2.21E-07, 2.96E-06, 1.11E-05, -2.77E-05, -0.000110969, 0.000240048, 0.000545759, -0.001557582, -0.001343505, 0.00950964, 0.017703103, 0.008650144, -0.00165215, -0.001379679, 0.00058102, 0.000199149, -0.000111584, -2.22E-05, 1.10E-05, 2.53E-06, -2.51E-07, -2.00E-07, -7.28E-08, -1.09E-08, 6.01E-09, 4.31E-09,
        2.27E-09, 4.64E-09, 5.63E-09, -1.44E-08, -8.15E-08, -2.16E-07, -2.02E-07, 3.18E-06, 1.11E-05, -3.06E-05, -0.00010998, 0.000260979, 0.000523636, -0.001644602, -0.001169811, 0.009936617, 0.017672083, 0.008219373, -0.001787326, -0.001289535, 0.000594324, 0.000179298, -0.000111263, -1.96E-05, 1.09E-05, 2.33E-06, -2.63E-07, -1.95E-07, -7.00E-08, -9.78E-09, 6.10E-09, 4.20E-09,
        2.29E-09, 4.76E-09, 5.46E-09, -1.57E-08, -8.46E-08, -2.21E-07, -1.80E-07, 3.42E-06, 1.11E-05, -3.36E-05, -0.000108502, 0.000282144, 0.000498421, -0.00172982, -0.000983105, 0.010360569, 0.017620465, 0.007789068, -0.001909942, -0.001199069, 0.000604864, 0.000159914, -0.000110559, -1.71E-05, 1.07E-05, 2.14E-06, -2.72E-07, -1.89E-07, -6.73E-08, -8.74E-09, 6.17E-09, 4.09E-09,
        2.30E-09, 4.87E-09, 5.28E-09, -1.70E-08, -8.77E-08, -2.27E-07, -1.56E-07, 3.67E-06, 1.10E-05, -3.68E-05, -0.000106507, 0.000303473, 0.000470051, -0.001812835, -0.000783356, 0.010780605, 0.017548376, 0.00736007, -0.002020186, -0.001108614, 0.000612739, 0.000141046, -0.000109498, -1.47E-05, 1.05E-05, 1.96E-06, -2.81E-07, -1.84E-07, -6.46E-08, -7.74E-09, 6.23E-09, 3.98E-09,
        2.32E-09, 4.98E-09, 5.07E-09, -1.84E-08, -9.08E-08, -2.32E-07, -1.29E-07, 3.92E-06, 1.08E-05, -3.99E-05, -0.000103972, 0.000324892, 0.000438473, -0.001893237, -0.000570569, 0.011195834, 0.017455989, 0.006933208, -0.002118275, -0.001018492, 0.000618054, 0.000122736, -0.000108107, -1.25E-05, 1.03E-05, 1.79E-06, -2.87E-07, -1.79E-07, -6.20E-08, -6.78E-09, 6.28E-09, 3.87E-09,
        2.35E-09, 5.09E-09, 4.84E-09, -1.98E-08, -9.40E-08, -2.38E-07, -9.81E-08, 4.19E-06, 1.07E-05, -4.32E-05, -0.000100872, 0.000346321, 0.000403644, -0.001970605, -0.000344782, 0.011605365, 0.017343528, 0.00650929, -0.002204453, -0.000929008, 0.000620916, 0.000105021, -0.000106412, -1.03E-05, 1.01E-05, 1.62E-06, -2.92E-07, -1.74E-07, -5.94E-08, -5.86E-09, 6.31E-09, 3.77E-09,
        2.38E-09, 5.20E-09, 4.59E-09, -2.13E-08, -9.73E-08, -2.43E-07, -6.39E-08, 4.46E-06, 1.04E-05, -4.66E-05, -9.72E-05, 0.000367675, 0.000365532, -0.002044512, -0.000106074, 0.01200831, 0.017211265, 0.006089103, -0.002278988, -0.00084045, 0.00062144, 8.79E-05, -0.000104439, -8.32E-06, 9.86E-06, 1.47E-06, -2.96E-07, -1.70E-07, -5.69E-08, -4.99E-09, 6.33E-09, 3.66E-09,
        2.42E-09, 5.30E-09, 4.31E-09, -2.28E-08, -1.01E-07, -2.49E-07, -2.58E-08, 4.74E-06, 1.01E-05, -5.00E-05, -9.29E-05, 0.000388868, 0.000324115, -0.002114521, 0.000145445, 0.012403787, 0.017059518, 0.005673413, -0.002342174, -0.00075309, 0.000619739, 7.15E-05, -0.000102213, -6.41E-06, 9.60E-06, 1.33E-06, -2.98E-07, -1.65E-07, -5.45E-08, -4.16E-09, 6.34E-09, 3.56E-09,
        2.46E-09, 5.41E-09, 4.01E-09, -2.44E-08, -1.04E-07, -2.54E-07, 1.64E-08, 5.02E-06, 9.72E-06, -5.34E-05, -8.79E-05, 0.000409805, 0.000279384, -0.002180191, 0.000409623, 0.012790924, 0.016888653, 0.00526296, -0.002394323, -0.000667186, 0.000615934, 5.58E-05, -9.98E-05, -4.62E-06, 9.33E-06, 1.19E-06, -3.00E-07, -1.60E-07, -5.21E-08, -3.37E-09, 6.33E-09, 3.46E-09,
        2.50E-09, 5.51E-09, 3.68E-09, -2.60E-08, -1.08E-07, -2.60E-07, 6.29E-08, 5.32E-06, 9.27E-06, -5.69E-05, -8.23E-05, 0.00043039, 0.000231341, -0.002241075, 0.000686271, 0.013168859, 0.016699079, 0.004858458, -0.002435771, -0.000582976, 0.000610144, 4.07E-05, -9.71E-05, -2.94E-06, 9.04E-06, 1.06E-06, -3.00E-07, -1.56E-07, -4.98E-08, -2.62E-09, 6.31E-09, 3.37E-09,
        2.55E-09, 5.61E-09, 3.33E-09, -2.77E-08, -1.11E-07, -2.65E-07, 1.14E-07, 5.62E-06, 8.74E-06, -6.05E-05, -7.61E-05, 0.000450523, 0.000180004, -0.002296725, 0.000975164, 0.013536744, 0.01649125, 0.004460593, -0.002466872, -0.000500683, 0.000602493, 2.64E-05, -9.43E-05, -1.37E-06, 8.74E-06, 9.40E-07, -3.00E-07, -1.51E-07, -4.75E-08, -1.90E-09, 6.29E-09, 3.28E-09,
        2.61E-09, 5.70E-09, 2.95E-09, -2.95E-08, -1.15E-07, -2.70E-07, 1.70E-07, 5.92E-06, 8.14E-06, -6.40E-05, -6.91E-05, 0.000470098, 0.000125402, -0.002346688, 0.001276037, 0.013893748, 0.016265663, 0.004070022, -0.002487996, -0.000420512, 0.000593102, 1.28E-05, -9.13E-05, 9.30E-08, 8.44E-06, 8.28E-07, -2.98E-07, -1.47E-07, -4.53E-08, -1.23E-09, 6.25E-09, 3.19E-09,
        2.67E-09, 5.79E-09, 2.54E-09, -3.12E-08, -1.19E-07, -2.75E-07, 2.31E-07, 6.23E-06, 7.46E-06, -6.76E-05, -6.15E-05, 0.000489008, 6.76E-05, -0.002390512, 0.001588586, 0.014239057, 0.016022853, 0.00368737, -0.00249953, -0.000342651, 0.000582097, -5.04E-08, -8.81E-05, 1.45E-06, 8.13E-06, 7.22E-07, -2.96E-07, -1.43E-07, -4.31E-08, -5.87E-10, 6.21E-09, 3.10E-09,
        2.73E-09, 5.87E-09, 2.10E-09, -3.31E-08, -1.22E-07, -2.79E-07, 2.98E-07, 6.55E-06, 6.69E-06, -7.11E-05, -5.31E-05, 0.000507142, 6.59E-06, -0.002427745, 0.001912472, 0.014571881, 0.015763397, 0.003313232, -0.002501874, -0.000267271, 0.0005696, -1.22E-05, -8.49E-05, 2.69E-06, 7.82E-06, 6.24E-07, -2.94E-07, -1.38E-07, -4.10E-08, 1.72E-11, 6.16E-09, 3.02E-09,
        2.80E-09, 5.95E-09, 1.63E-09, -3.50E-08, -1.26E-07, -2.83E-07, 3.70E-07, 6.86E-06, 5.83E-06, -7.46E-05, -4.40E-05, 0.000524386, -5.75E-05, -0.00245794, 0.002247313, 0.014891453, 0.015487908, 0.002948165, -0.00249544, -0.000194528, 0.000555735, -2.35E-05, -8.15E-05, 3.84E-06, 7.50E-06, 5.33E-07, -2.91E-07, -1.34E-07, -3.89E-08, 5.87E-10, 6.09E-09, 2.94E-09,
        2.87E-09, 6.03E-09, 1.12E-09, -3.69E-08, -1.30E-07, -2.87E-07, 4.48E-07, 7.18E-06, 4.88E-06, -7.81E-05, -3.41E-05, 0.000540623, -0.000124558, -0.002480651, 0.002592696, 0.015197033, 0.015197033, 0.002592696, -0.002480651, -0.000124558, 0.000540623, -3.41E-05, -7.81E-05, 4.88E-06, 7.18E-06, 4.48E-07, -2.87E-07, -1.30E-07, -3.69E-08, 1.12E-09, 6.03E-09, 2.87E-09,
        2.94E-09, 6.09E-09, 5.87E-10, -3.89E-08, -1.34E-07, -2.91E-07, 5.33E-07, 7.50E-06, 3.84E-06, -8.15E-05, -2.35E-05, 0.000555735, -0.000194528, -0.00249544, 0.002948165, 0.015487908, 0.014891453, 0.002247313, -0.00245794, -5.75E-05, 0.000524386, -4.40E-05, -7.46E-05, 5.83E-06, 6.86E-06, 3.70E-07, -2.83E-07, -1.26E-07, -3.50E-08, 1.63E-09, 5.95E-09, 2.80E-09,
        3.02E-09, 6.16E-09, 1.72E-11, -4.10E-08, -1.38E-07, -2.94E-07, 6.24E-07, 7.82E-06, 2.69E-06, -8.49E-05, -1.22E-05, 0.0005696, -0.000267271, -0.002501874, 0.003313232, 0.015763397, 0.014571881, 0.001912472, -0.002427745, 6.59E-06, 0.000507142, -5.31E-05, -7.11E-05, 6.69E-06, 6.55E-06, 2.98E-07, -2.79E-07, -1.22E-07, -3.31E-08, 2.10E-09, 5.87E-09, 2.73E-09,
        3.10E-09, 6.21E-09, -5.87E-10, -4.31E-08, -1.43E-07, -2.96E-07, 7.22E-07, 8.13E-06, 1.45E-06, -8.81E-05, -5.04E-08, 0.000582097, -0.000342651, -0.00249953, 0.00368737, 0.016022853, 0.014239057, 0.001588586, -0.002390512, 6.76E-05, 0.000489008, -6.15E-05, -6.76E-05, 7.46E-06, 6.23E-06, 2.31E-07, -2.75E-07, -1.19E-07, -3.12E-08, 2.54E-09, 5.79E-09, 2.67E-09,
        3.19E-09, 6.25E-09, -1.23E-09, -4.53E-08, -1.47E-07, -2.98E-07, 8.28E-07, 8.44E-06, 9.30E-08, -9.13E-05, 1.28E-05, 0.000593102, -0.000420512, -0.002487996, 0.004070022, 0.016265663, 0.013893748, 0.001276037, -0.002346688, 0.000125402, 0.000470098, -6.91E-05, -6.40E-05, 8.14E-06, 5.92E-06, 1.70E-07, -2.70E-07, -1.15E-07, -2.95E-08, 2.95E-09, 5.70E-09, 2.61E-09,
        3.28E-09, 6.29E-09, -1.90E-09, -4.75E-08, -1.51E-07, -3.00E-07, 9.40E-07, 8.74E-06, -1.37E-06, -9.43E-05, 2.64E-05, 0.000602493, -0.000500683, -0.002466872, 0.004460593, 0.01649125, 0.013536744, 0.000975164, -0.002296725, 0.000180004, 0.000450523, -7.61E-05, -6.05E-05, 8.74E-06, 5.62E-06, 1.14E-07, -2.65E-07, -1.11E-07, -2.77E-08, 3.33E-09, 5.61E-09, 2.55E-09,
        3.37E-09, 6.31E-09, -2.62E-09, -4.98E-08, -1.56E-07, -3.00E-07, 1.06E-06, 9.04E-06, -2.94E-06, -9.71E-05, 4.07E-05, 0.000610144, -0.000582976, -0.002435771, 0.004858458, 0.016699079, 0.013168859, 0.000686271, -0.002241075, 0.000231341, 0.00043039, -8.23E-05, -5.69E-05, 9.27E-06, 5.32E-06, 6.29E-08, -2.60E-07, -1.08E-07, -2.60E-08, 3.68E-09, 5.51E-09, 2.50E-09,
        3.46E-09, 6.33E-09, -3.37E-09, -5.21E-08, -1.60E-07, -3.00E-07, 1.19E-06, 9.33E-06, -4.62E-06, -9.98E-05, 5.58E-05, 0.000615934, -0.000667186, -0.002394323, 0.00526296, 0.016888653, 0.012790924, 0.000409623, -0.002180191, 0.000279384, 0.000409805, -8.79E-05, -5.34E-05, 9.72E-06, 5.02E-06, 1.64E-08, -2.54E-07, -1.04E-07, -2.44E-08, 4.01E-09, 5.41E-09, 2.46E-09,
        3.56E-09, 6.34E-09, -4.16E-09, -5.45E-08, -1.65E-07, -2.98E-07, 1.33E-06, 9.60E-06, -6.41E-06, -0.000102213, 7.15E-05, 0.000619739, -0.00075309, -0.002342174, 0.005673413, 0.017059518, 0.012403787, 0.000145445, -0.002114521, 0.000324115, 0.000388868, -9.29E-05, -5.00E-05, 1.01E-05, 4.74E-06, -2.58E-08, -2.49E-07, -1.01E-07, -2.28E-08, 4.31E-09, 5.30E-09, 2.42E-09,
        3.66E-09, 6.33E-09, -4.99E-09, -5.69E-08, -1.70E-07, -2.96E-07, 1.47E-06, 9.86E-06, -8.32E-06, -0.000104439, 8.79E-05, 0.00062144, -0.00084045, -0.002278988, 0.006089103, 0.017211265, 0.01200831, -0.000106074, -0.002044512, 0.000365532, 0.000367675, -9.72E-05, -4.66E-05, 1.04E-05, 4.46E-06, -6.39E-08, -2.43E-07, -9.73E-08, -2.13E-08, 4.59E-09, 5.20E-09, 2.38E-09,
        3.77E-09, 6.31E-09, -5.86E-09, -5.94E-08, -1.74E-07, -2.92E-07, 1.62E-06, 1.01E-05, -1.03E-05, -0.000106412, 0.000105021, 0.000620916, -0.000929008, -0.002204453, 0.00650929, 0.017343528, 0.011605365, -0.000344782, -0.001970605, 0.000403644, 0.000346321, -0.000100872, -4.32E-05, 1.07E-05, 4.19E-06, -9.81E-08, -2.38E-07, -9.40E-08, -1.98E-08, 4.84E-09, 5.09E-09, 2.35E-09,
        3.87E-09, 6.28E-09, -6.78E-09, -6.20E-08, -1.79E-07, -2.87E-07, 1.79E-06, 1.03E-05, -1.25E-05, -0.000108107, 0.000122736, 0.000618054, -0.001018492, -0.002118275, 0.006933208, 0.017455989, 0.011195834, -0.000570569, -0.001893237, 0.000438473, 0.000324892, -0.000103972, -3.99E-05, 1.08E-05, 3.92E-06, -1.29E-07, -2.32E-07, -9.08E-08, -1.84E-08, 5.07E-09, 4.98E-09, 2.32E-09,
        3.98E-09, 6.23E-09, -7.74E-09, -6.46E-08, -1.84E-07, -2.81E-07, 1.96E-06, 1.05E-05, -1.47E-05, -0.000109498, 0.000141046, 0.000612739, -0.001108614, -0.002020186, 0.00736007, 0.017548376, 0.010780605, -0.000783356, -0.001812835, 0.000470051, 0.000303473, -0.000106507, -3.68E-05, 1.10E-05, 3.67E-06, -1.56E-07, -2.27E-07, -8.77E-08, -1.70E-08, 5.28E-09, 4.87E-09, 2.30E-09,
        4.09E-09, 6.17E-09, -8.74E-09, -6.73E-08, -1.89E-07, -2.72E-07, 2.14E-06, 1.07E-05, -1.71E-05, -0.000110559, 0.000159914, 0.000604864, -0.001199069, -0.001909942, 0.007789068, 0.017620465, 0.010360569, -0.000983105, -0.00172982, 0.000498421, 0.000282144, -0.000108502, -3.36E-05, 1.11E-05, 3.42E-06, -1.80E-07, -2.21E-07, -8.46E-08, -1.57E-08, 5.46E-09, 4.76E-09, 2.29E-09,
        4.20E-09, 6.10E-09, -9.78E-09, -7.00E-08, -1.95E-07, -2.63E-07, 2.33E-06, 1.09E-05, -1.96E-05, -0.000111263, 0.000179298, 0.000594324, -0.001289535, -0.001787326, 0.008219373, 0.017672083, 0.009936617, -0.001169811, -0.001644602, 0.000523636, 0.000260979, -0.00010998, -3.06E-05, 1.11E-05, 3.18E-06, -2.02E-07, -2.16E-07, -8.15E-08, -1.44E-08, 5.63E-09, 4.64E-09, 2.27E-09,
        4.31E-09, 6.01E-09, -1.09E-08, -7.28E-08, -2.00E-07, -2.51E-07, 2.53E-06, 1.10E-05, -2.22E-05, -0.000111584, 0.000199149, 0.00058102, -0.001379679, -0.00165215, 0.008650144, 0.017703103, 0.00950964, -0.001343505, -0.001557582, 0.000545759, 0.000240048, -0.000110969, -2.77E-05, 1.11E-05, 2.96E-06, -2.21E-07, -2.11E-07, -7.85E-08, -1.32E-08, 5.77E-09, 4.53E-09, 2.27E-09];

    foreach(a, b; lockstep(coeff.flattened, testResult.sliced(coeff.shape).flattened))
        assert(approxEqual(a, b));

    size_t delay = 64 * (32 - 1) / 2;
    size_t len = 32*32*2;
    auto chirp = sequence!"(a[0]*n*(n-1)/2) % 1"(1.0/len).unaryFun!(a => std.algorithm.map!(a => std.complex.expi(-a * 2 * PI))(a).take(len).array());
    //auto chirp = sequence!"a[0]*n"(0.125).map!(a => std.complex.expi(a * 2 * PI)).take(len).array();
    auto analysed = analysis(coeff, chirp);

    //File file = File("foo.csv","w");
    //file.writefln("%(%(%s,%),\n%)", analysed.mapSlice!(a => 10*log10(a.re^^2 + a.im^^2)).mapSlice!(a => a < -50 ? -50 : a));
    
    auto synthesized = synthesis(coeff, analysed)[delay .. $];
    chirp = chirp[0 .. synthesized.length];

    real sum = 0, P = 0;
    foreach(e; synthesized.zip(chirp)){
        assert(approxEqual(e[0].re, e[1].re));
        assert(approxEqual(e[0].im, e[1].im));
        sum += (e[0] - e[1]).sqAbs;
        P += e[0].sqAbs;
    }

    assert(10*log10(sum / P) < -130);
}


private
auto matlabReshape(S)(S s, size_t n, size_t m)
{
    int e;
    auto rs = s.reshape([m, n], e).universal.transposed;
    return rs;
}
