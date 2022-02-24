module dffdd.detector.gabp;

import std.algorithm;
import std.complex;
import std.math;
import std.traits;

import mir.ndslice;
import mir.blas;

import dffdd.detector.primitives;
import dffdd.mod;
import dffdd.math;


struct GaBPWorkspace(C, Mod)
if(is(Mod == BPSK!C) || is(Mod == QPSK!C))
{
    this(Mod mod, size_t ntxvec, size_t nrxvec)
    {
        _mod = mod;
        _ntxvec = ntxvec;
        _nrxvec = nrxvec;
        _chMat = slice!C(nrxvec, ntxvec);
        _chMatSQ = slice!R(nrxvec, ntxvec);
        _recvY = new C[nrxvec];
        _sigma2 = new R[nrxvec];
        _lastLLR = new SoftSymType[ntxvec];
        _softSym = new SoftSymType[ntxvec];
        _expSoftSymPower = new R[ntxvec];

        this.clear();
    }


    void channelMatrix(Matrix)(in Matrix channel)
    if(is(typeof(channel[0][0]) : C))
    {
        foreach(i; 0 .. _nrxvec)
            foreach(j; 0 .. _ntxvec) {
                immutable h = channel[i][j];
                _chMat[i][j] = h;
                _chMatSQ[i][j] = h.sqAbs;
            }
    }


    void clear()
    {
        _recvY[] = C(0);
        _sigma2[] = R(0);
        _lastLLR[] = SoftSymType(0);
        _softSym[] = SoftSymType(0);
        _expSoftSymPower[] = R(1);
    }


    ref Bit[] decide(return ref Bit[] bits)
    {
        static if(is(SoftSymType == R))
        {
            bits.length = _ntxvec;
            foreach(i; 0 .. _ntxvec) {
                bits[i] = _softSym[i] > 0 ? 0 : 1;
            }
        }
        else
        {
            bits.length = _ntxvec * 2;
            foreach(i; 0 .. _ntxvec) {
                bits[i*2 + 0] = _softSym[i].re > 0 ? 0 : 1;
                bits[i*2 + 1] = _softSym[i].im > 0 ? 0 : 1;
            }
        }

        return bits;
    }


    inout(R)[] lastLLR() inout
    {
        static if(is(SoftSymType == R))
            enum BaseElemSize = 1;
        else
            enum BaseElemSize = 2;

        return (cast(inout(R)*)_lastLLR.ptr)[0 .. _lastLLR.length * BaseElemSize];
    }


    const(SoftSymType)[] calcSoftSymbol()
    {
        foreach(i; 0 .. _ntxvec) {
            static if(is(SoftSymType == R))
            {
                _softSym[i] = tanh(_lastLLR[i]/2);
                _expSoftSymPower[i] = 1 - _softSym[i]^^2;
            }
            else
            {
                _softSym[i] = C(tanh(_lastLLR[i].re/2), tanh(_lastLLR[i].im/2)) * SQRT1_2;
                _expSoftSymPower[i] = 1 - _softSym[i].sqAbs;
            }
        }

        return _softSym;
    }


    // void lastLLR(in R[] llrs)
    // {
    //     _lastLLR[] = llrs[];
    // }


  private:
    alias R = typeof(C.init.re);
    alias SoftSymType = Select!(is(Mod == BPSK!C), R, C);

    Mod _mod;

    size_t _ntxvec, _nrxvec;

    Slice!(C*, 2, Contiguous) _chMat;
    Slice!(R*, 2, Contiguous) _chMatSQ;
    Slice!(C*, 2, Contiguous) _chMatQ, _chMatR;
    C[] _recvY;
    R[] _sigma2;
    SoftSymType[] _lastLLR;
    SoftSymType[] _softSym;
    R[] _expSoftSymPower;
}



void iterateGaBPDetection(C, Mod)(ref GaBPWorkspace!(C, Mod) ws, in C[] recvY, double N0, double dummping, double LLR_MAX)
if(is(Mod == BPSK!C) || is(Mod == QPSK!C))
{
    alias SSType = ws.SoftSymType;
    alias R = typeof(C.init.re);

    ws._recvY[] = recvY[];
    ws._sigma2[] = N0;

    ws.calcSoftSymbol();

    mir.blas.gemv(C(-1), ws._chMat, ws._softSym.sliced, C(1), ws._recvY.sliced);
    mir.blas.gemv(R(1), ws._chMatSQ, ws._expSoftSymPower.sliced, R(1), ws._sigma2.sliced);
    // foreach(j; 0 .. ws._nrxvec) {
    //     foreach(i; 0 .. ws._ntxvec) {
    //         ws._recvY[j] -= ws._chMat[j, i] * ws._softSym[i];
    //         ws._sigma2[j] += ws._chMatSQ[j, i] * ws._expSoftSymPower[i];
    //     }
    // }


    foreach(i; 0 .. ws._ntxvec) {
        SSType newLLR = SSType(0);
        foreach(j; 0 .. ws._nrxvec) {
            immutable y_ij = ws._recvY[j] + ws._chMat[j, i] * ws._softSym[i];
            immutable sig_ij = ws._sigma2[j] - ws._chMatSQ[j, i] * ws._expSoftSymPower[i];

            static if(is(SSType == R))
                newLLR += 4 * (y_ij * ws._chMat[j, i].conj).re / sig_ij;
            else
                newLLR += 4 * (y_ij * ws._chMat[j, i].conj) / sig_ij;
        }

        static if(is(SSType == R))
            newLLR = min(max(newLLR, -LLR_MAX), +LLR_MAX);
        else
        {
            newLLR.re = min(max(newLLR.re, -LLR_MAX), +LLR_MAX);
            newLLR.im = min(max(newLLR.im, -LLR_MAX), +LLR_MAX);
        }

        newLLR = newLLR * (1 - dummping) + ws._lastLLR[i] * dummping;
        ws._lastLLR[i] = newLLR;
    }

}


enum GaBPDetectorMode
{
    toBit,          // C -> Bit
    toSoftSymbol,   // C -> C
    toLLR,          // C -> R
    toP0P1          // C -> R
}


class GaBPDetector(C, Mod, GaBPDetectorMode mode)
    : IDetector!(C, Select!(mode == GaBPDetectorMode.toBit, Bit,
                    Select!(mode == GaBPDetectorMode.toSoftSymbol, C, typeof(C.init.re))))
{
    alias InputElementType = C;

  static if(mode == GaBPDetectorMode.toBit)
    alias OutputElementType = Bit;
  else static if(mode == GaBPDetectorMode.toSoftSymbol)
    alias OutputElementType = C;
  else
    alias OutputElementType = typeof(C.init.re);



    this(Matrix)(Mod mod, size_t ntxvec, size_t nrxvec, Matrix chMat, double N0, size_t maxIter = 20, double dummping = 0.75, double LLR_MAX = 4)
    if(is(typeof(chMat[0][0]) : C))
    {
        _ws = GaBPWorkspace!(C, Mod)(mod, ntxvec, nrxvec);
        _ws.channelMatrix = chMat;
        _ntxvec = ntxvec;
        _nrxvec = nrxvec;
        _maxIter = maxIter;
        _N0 = N0;
        _dummping = dummping;
        _LLR_MAX = LLR_MAX;
    }


    size_t inputLength() const
    {
        return _nrxvec;
    }


    size_t outputLength() const
    {
        return _ntxvec;
    }


    OutputElementType[] detect(in C[] received, return ref OutputElementType[] detected)
    in(received.length % _nrxvec == 0)
    do {
        R[] null_;
        this.detectImpl(received, null_, detected);
        return detected;
    }


    E[] detect(E)(in C[] received, return ref E[] detected)
    if(is(E : OutputElementType) && !is(E == OutputElementType))
    in(received.length % _nrxvec == 0)
    do {
        R[] null_;
        this.detectImpl(received, null_, detected);
        return detected;
    }


    E2[] detect(E1, E2)(in C[] received, in E1[] priorLLR, return ref E2[] detected)
    in(received.length % _nrxvec == 0)
    in(priorLLR.length == received.length / _nrxvec * _ntxvec * 2)
    do {
        this.detectImpl(received, priorLLR, detected);
        return detected;
    }


  private:
    GaBPWorkspace!(C, Mod) _ws;
    size_t _ntxvec;
    size_t _nrxvec;
    size_t _maxIter;
    double _N0;
    double _dummping;
    double _LLR_MAX;


    alias R = typeof(C.init.re);


    enum BaseElemSize = is(Mod == BPSK!C) ? 1 : 2;


    void[] detectImpl(E1, E2)(in C[] received, in E1[] priorLLR, ref E2[] detected)
    if(is(E1 : R) && is(E2 : OutputElementType))
    {
        immutable nblock = received.length / _nrxvec;

        static if(mode == GaBPDetectorMode.toSoftSymbol)
            detected.length = nblock * _ntxvec;
        else    // if Bit, LLR, or P0P1
            detected.length = nblock * _ntxvec * BaseElemSize;

        immutable blockSize = detected.length / nblock;

        foreach(i; 0 .. nblock) {
            if(priorLLR is null)
                this.detectEach(received[_nrxvec * i .. _nrxvec * (i+1)], cast(R[])null, detected[blockSize * i .. blockSize * (i+1)]);
            else
                this.detectEach(received[_nrxvec * i .. _nrxvec * (i+1)], priorLLR[_ntxvec * i * 2 .. _ntxvec * (i+1) * 2], detected[blockSize * i .. blockSize * (i+1)]);
        }

        return detected;
    }


    void detectEachImpl(E1)(in C[] received, in E1[] priorLLR)
    if(is(E1 : R))
    {
        _ws.clear();
        if(priorLLR !is null) {
            foreach(i; 0 .. _ntxvec * 2)
                _ws.lastLLR[i] = priorLLR[i];
        }

        foreach(i; 0 .. _maxIter)
            iterateGaBPDetection(_ws, received, _N0, _dummping, _LLR_MAX);
    }


  static if(mode == GaBPDetectorMode.toBit)
  {
    void detectEach(E1, E2)(in C[] received, in E1[] priorLLR, E2[] detected)
    if(is(E1 : R) && is(E2 : Bit))
    {
        static Bit[] tmp; 

        this.detectEachImpl(received, priorLLR);
        _ws.decide(tmp);

        foreach(i; 0 .. _ntxvec * 2)
            detected[i] = tmp[i];
    }
  }
  else static if(mode == GaBPDetectorMode.toSoftSymbol)
  {
    void detectEach(E1, E2)(in C[] received, in E1[] priorLLR, E2[] detected)
    if(is(E1 : R) && is(E2 : C))
    {
        this.detectEachImpl(received, priorLLR);
        auto ss = _ws.calcSoftSymbol();

        foreach(i; 0 .. _ntxvec)
            detected[i] = ss[i];
    }
  }
  else
  {
    void detectEach(E1, E2)(in C[] received, in E1[] priorLLR, E2[] detected)
    if(is(E1 : R) && is(E2 : R))
    {
        this.detectEachImpl(received, priorLLR);
        auto llr = _ws.lastLLR;

        static if(mode == GaBPDetectorMode.toLLR)
        {
            foreach(i; 0 .. _ntxvec * 2)
                detected[i] = llr[i];
        }
        else static if(mode == GaBPDetectorMode.toP0P1)
        {
            foreach(i; 0 .. _ntxvec * 2)
                detected[i] = exp(llr[i]);
        }
        else static assert(0, "Unsupported mode: %s".format(mode));
    }
  }
}


// unittest
// {
//     import dffdd.mod.bpsk;
//     alias C = Complex!double;

//     auto mod = BPSK!C();
//     auto gabp = new GaBPDetector!(C, typeof(mod), GaBPDetectorMode.toP0P1)(mod, 1, 1, [[C(1)]], 0);
// }


// struct GaussianBPDetector
// {
//     this(M, Mod)(size_t numInput, size_t numOutput, M channelMatrix, Mod mod)
//     {
        
//     }
// }
