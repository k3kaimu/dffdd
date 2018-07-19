module preamble;

import std.algorithm : max, min;
import std.complex;
import std.typecons;
import std.stdio;
import std.math;

import blas;
import spec;


struct Preambles
{
  shared static
  {
    // shared padding
    immutable(Complex!float[]) shdPAD;

    // 時間ドメイン(1023 * nOSサンプル)
    immutable(Complex!float[]) fstEST;
    immutable(Complex!float[]) fstRTS;
    immutable(Complex!float[]) fstCTS;
    immutable(Complex!float[]) fstDAT;
    immutable(Complex!float[]) fstACK;

    immutable(Complex!float[]) sndEST;
    immutable(Complex!float[]) sndRTS;
    immutable(Complex!float[]) sndCTS;
    immutable(Complex!float[]) sndDAT;
    immutable(Complex!float[]) sndACK;

    // 周波数ドメイン(2048 * nOSサンプル = FFT(code(1023) + zero(1025)) * nOS)
    immutable(Complex!float[]) fstESTFreq;
    immutable(Complex!float[]) fstRTSFreq;
    immutable(Complex!float[]) fstCTSFreq;
    immutable(Complex!float[]) fstDATFreq;
    immutable(Complex!float[]) fstACKFreq;

    immutable(Complex!float[]) sndESTFreq;
    immutable(Complex!float[]) sndRTSFreq;
    immutable(Complex!float[]) sndCTSFreq;
    immutable(Complex!float[]) sndDATFreq;
    immutable(Complex!float[]) sndACKFreq;
  }

    enum uint prnOfSHDPAD = 100;

    enum uint prnOfFSTEST = 1;
    enum uint prnOfFSTRTS = 2;
    enum uint prnOfFSTCTS = 3;
    enum uint prnOfFSTDAT = 4;
    enum uint prnOfFSTACK = 5;

    enum uint prnOfSNDEST = 11;
    enum uint prnOfSNDRTS = 12;
    enum uint prnOfSNDCTS = 13;
    enum uint prnOfSNDDAT = 14;
    enum uint prnOfSNDACK = 15;

    shared static this()
    {
        import std.range : take;
        import dffdd.gps.code : L1CACode;

        immutable(Complex!float[]) makeList(uint prn)
        {
            Complex!float[] ret;
            auto code = L1CACode(prn);

            foreach(byte e; code.take(1023)){
                foreach(i; 0 .. Constant.nOverSampling)
                    ret ~= Complex!float(e, 0) * 0.1f;
            }

            return cast(immutable)ret;
        }

        Preambles.shdPAD = makeList(Preambles.prnOfSHDPAD);

        Preambles.fstEST = makeList(Preambles.prnOfFSTEST);
        Preambles.fstRTS = makeList(Preambles.prnOfFSTRTS);
        Preambles.fstCTS = makeList(Preambles.prnOfFSTCTS);
        Preambles.fstDAT = makeList(Preambles.prnOfFSTDAT);
        Preambles.fstACK = makeList(Preambles.prnOfFSTACK);

        Preambles.sndEST = makeList(Preambles.prnOfSNDEST);
        Preambles.sndRTS = makeList(Preambles.prnOfSNDRTS);
        Preambles.sndCTS = makeList(Preambles.prnOfSNDCTS);
        Preambles.sndDAT = makeList(Preambles.prnOfSNDDAT);
        Preambles.sndACK = makeList(Preambles.prnOfSNDACK);

        immutable(Complex!float[]) makeFreqDomain(immutable(Complex!float)[] buffer)
        {
            import std.numeric : fft;

            foreach(i; 0 .. 1025 * Constant.nOverSampling)
                buffer ~= Complex!float(0, 0);

            return fft!float(buffer).dup;
        }

        Preambles.fstESTFreq = makeFreqDomain(Preambles.fstEST);
        Preambles.fstRTSFreq = makeFreqDomain(Preambles.fstRTS);
        Preambles.fstCTSFreq = makeFreqDomain(Preambles.fstCTS);
        Preambles.fstDATFreq = makeFreqDomain(Preambles.fstDAT);
        Preambles.fstACKFreq = makeFreqDomain(Preambles.fstACK);

        Preambles.sndESTFreq = makeFreqDomain(Preambles.sndEST);
        Preambles.sndRTSFreq = makeFreqDomain(Preambles.sndRTS);
        Preambles.sndCTSFreq = makeFreqDomain(Preambles.sndCTS);
        Preambles.sndDATFreq = makeFreqDomain(Preambles.sndDAT);
        Preambles.sndACKFreq = makeFreqDomain(Preambles.sndACK);
    }


    this(bool isSecond)
    {
        if(!isSecond)
        {
            myEST = fstEST;
            myRTS = fstRTS;
            myCTS = fstCTS;
            myDAT = fstDAT;
            myACK = fstACK;

            tgEST = sndEST;
            tgRTS = sndRTS;
            tgCTS = sndCTS;
            tgDAT = sndDAT;
            tgACK = sndACK;

            myESTFreq = fstESTFreq;
            myRTSFreq = fstRTSFreq;
            myCTSFreq = fstCTSFreq;
            myDATFreq = fstDATFreq;
            myACKFreq = fstACKFreq;

            tgESTFreq = sndESTFreq;
            tgRTSFreq = sndRTSFreq;
            tgCTSFreq = sndCTSFreq;
            tgDATFreq = sndDATFreq;
            tgACKFreq = sndACKFreq;
        }
        else
        {
            myEST = sndEST;
            myRTS = sndRTS;
            myCTS = sndCTS;
            myDAT = sndDAT;
            myACK = sndACK;

            tgEST = fstEST;
            tgRTS = fstRTS;
            tgCTS = fstCTS;
            tgDAT = fstDAT;
            tgACK = fstACK;

            myESTFreq = sndESTFreq;
            myRTSFreq = sndRTSFreq;
            myCTSFreq = sndCTSFreq;
            myDATFreq = sndDATFreq;
            myACKFreq = sndACKFreq;

            tgESTFreq = fstESTFreq;
            tgRTSFreq = fstRTSFreq;
            tgCTSFreq = fstCTSFreq;
            tgDATFreq = fstDATFreq;
            tgACKFreq = fstACKFreq;
        }
    }


    immutable(Complex!float)[] myEST;
    immutable(Complex!float)[] myRTS;
    immutable(Complex!float)[] myCTS;
    immutable(Complex!float)[] myDAT;
    immutable(Complex!float)[] myACK;

    immutable(Complex!float)[] tgEST;
    immutable(Complex!float)[] tgRTS;
    immutable(Complex!float)[] tgCTS;
    immutable(Complex!float)[] tgDAT;
    immutable(Complex!float)[] tgACK;

    // 周波数ドメイン(2048 * nOSサンプル = FFT(code(1023) + zero(1025)) * nOS)
    immutable(Complex!float)[] myESTFreq;
    immutable(Complex!float)[] myRTSFreq;
    immutable(Complex!float)[] myCTSFreq;
    immutable(Complex!float)[] myDATFreq;
    immutable(Complex!float)[] myACKFreq;

    immutable(Complex!float)[] tgESTFreq;
    immutable(Complex!float)[] tgRTSFreq;
    immutable(Complex!float)[] tgCTSFreq;
    immutable(Complex!float)[] tgDATFreq;
    immutable(Complex!float)[] tgACKFreq;
}


final class PreambleDetector
{
    import dffdd.utils.fft;

    this(Preambles p)
    {
        this(p, 0);
    }


    this(Preambles p, real offsetFreq)
    {
        this._p = p;
        this._fftObj = FFTWObject!(Complex!float)(1024 * Constant.nOverSampling * 2);
        // this._offsetFreq = offsetFreq;
        this._reSampling();
    }


    // void offsetFreq(real freq) @property 
    // {
    //     // if(_offsetFreq != freq){
    //     //     _offsetFreq = freq;
    //     //     _reSampling();
    //     // }else
    //     //     _offsetFreq = freq;
    // }


    // real offsetFreq() const @property { return _offsetFreq; }


    Tuple!(bool, size_t) detect(string tgt)(in Complex!float[] haystack)
    if(tgt.length == 5
        && (tgt[0 .. 2] == "My" || tgt[0 .. 2] == "Tg")
        && (tgt[2 .. 5] == "EST" || tgt[2 .. 5] == "RTS" || tgt[2 .. 5] == "CTS"
            || tgt[2 .. 5] == "DAT" || tgt[2 .. 5] == "ACK"
    ))
    in{
        assert(haystack.length % (Constant.nOverSampling * 1024) == 0);
        assert(haystack.length >= Constant.nOverSampling * 1024 * 2);
    }
    body{
        immutable size_t sizeOfFFT = Constant.nOverSampling * 1024 * 2;
        immutable(Complex!float)[] needle;

        static if(tgt[0 .. 2] == "My")
            needle = mixin(`_p.my` ~ tgt[2 .. 5] ~ `Freq`);
        else
        {
            switch(tgt[2 .. 5]){
                case "EST": needle = _tgPreambleFreq[0]; break;
                case "RTS": needle = _tgPreambleFreq[1]; break;
                case "CTS": needle = _tgPreambleFreq[2]; break;
                case "DAT": needle = _tgPreambleFreq[3]; break;
                case "ACK": needle = _tgPreambleFreq[4]; break;
                default: assert(0);
            }
        }

        //_fftObj.inputs!float[] = haystack[];
        size_t maxTryCount = (haystack.length / (sizeOfFFT/2)) - 1;
        foreach(tryIndex; 0 .. maxTryCount){
            immutable baseIdx = tryIndex * (sizeOfFFT/2);

            auto ips = _fftObj.inputs!float;
            auto ops = _fftObj.outputs!float;

            ips[] = haystack[baseIdx .. baseIdx + sizeOfFFT];
            _fftObj.fft!float();
            ips[] = ops[];

            foreach(i; 0 .. sizeOfFFT)
                ips[i] *= -needle[i].conj;

            _fftObj.ifft!float();
            ptrdiff_t idxp1 = BLAS.ixamax(ops[0 .. 1023 * Constant.nOverSampling]);
            assert(idxp1 >= 0);

            immutable maxP1 = ops[idxp1].sqAbs;
            immutable size_t exclude = Constant.nOverSampling*4;

            ptrdiff_t idxp2;
            real maxP2;

            if(idxp1 < exclude || (1023 - idxp1) < exclude){
                ptrdiff_t sidx, eidx;
                if(idxp1 < exclude){
                    sidx = idxp1 + exclude;
                    eidx = 1023 - exclude + idxp1;
                }else{
                    sidx = exclude - (1023 - idxp1);
                    eidx = idxp1 - exclude;
                }

                idxp2 = BLAS.ixamax(ops[sidx .. eidx]);
                maxP2 = ops[idxp2].sqAbs;
            }else{
                immutable exidx1 = idxp1 - cast(ptrdiff_t)exclude;
                immutable exidx2 = idxp1 + exclude;
                ptrdiff_t idxp21 = BLAS.ixamax(ops[0 .. exidx1]);
                ptrdiff_t idxp22 = BLAS.ixamax(ops[exidx2 .. $]);
                assert(idxp21 >= 0);
                assert(idxp22 >= 0);

                immutable maxP21 = ops[idxp21].sqAbs,
                          maxP22 = ops[idxp22 + exidx2].sqAbs;

                if(maxP21 > maxP22){
                    idxp2 = idxp21;
                    maxP2 = maxP21;
                }else{
                    idxp2 = idxp22;
                    maxP2 = maxP22;
                }
            }

             if(maxP1 > 1){
                  writefln("%s : %s", idxp1, idxp2);
                  writeln(maxP1, " / ", maxP2, " = ", maxP1 / maxP2);
             }

            if(maxP1 / maxP2 > Constant.PreambleDetector.acqTH){
                // writeln(maxP1, " / ", maxP2, " = ", maxP1 / maxP2);
                return typeof(return)(true, baseIdx + idxp1);
            }
        }

        //writeln("Failed");
        return typeof(return)(false, 0);
    }


  private:
    Preambles _p;
    FFTWObject!(Complex!float) _fftObj;
    // real _offsetFreq = 0;
    immutable(Complex!float)[][5] _tgPreambleFreq;   // [est, rts, cts, dat, ack]


    void _reSampling()
    {
        foreach(codeIdx, s; [_p.tgEST, _p.tgRTS, _p.tgCTS, _p.tgDAT, _p.tgACK]){
            import carbon.stream;

            auto ips = _fftObj.inputs!float;
            ips[0 .. 1023 * Constant.nOverSampling] = s[];

            // zero-padding
            ips[1023 * Constant.nOverSampling .. $] = Complex!float(0, 0);

            // mul nco
            // auto nco = lutNCO!(std.complex.expi, 1024*1024)(_offsetFreq, 1.0L / Constant.samplingFreq);
            // nco.readOp!"*"(ips[0 .. 1023 * Constant.nOverSampling]);

            // fft
            _fftObj.fft!float();
            _tgPreambleFreq[codeIdx] = _fftObj.outputs!float().dup;
        }
    }
}


unittest
{
    immutable testIndex = [0, 1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400, 1023 * Constant.nOverSampling - 9, 1023 * Constant.nOverSampling - 8, 1023 * Constant.nOverSampling - 7, 1023 * Constant.nOverSampling - 6, 1023 * Constant.nOverSampling - 5, 1023 * Constant.nOverSampling - 4, 1023 * Constant.nOverSampling - 3, 1023 * Constant.nOverSampling - 2, 1023 * Constant.nOverSampling - 1];
    immutable testFreq = [0];

    Preambles preambles = Preambles(false);
    //import std.datetime;
    foreach(freq; testFreq){
        auto detector = new PreambleDetector(preambles, freq);

        auto copied = preambles.tgRTS.dup;
        foreach(ref e; copied) e *= Complex!float(1, 1);

        //StopWatch sw;
        //sw.start();
        foreach(offsetIdx; testIndex){
            foreach(offsetIdx2; testIndex){
                foreach(b; [true, false]){
                    Complex!float[] recv = new Complex!float[2048 * Constant.nOverSampling];
                    recv[] = Complex!float(0, 0);

                    if(b) foreach(i; 0 .. 1023*Constant.nOverSampling)
                        recv[offsetIdx + i] = copied[i];

                    foreach(i; 0 .. 1023*Constant.nOverSampling)
                        recv[offsetIdx2 + i] += preambles.myRTS[i] * Complex!float(1, 1);

                    auto res = detector.detect!"TgRTS"(recv);
                    assert(res[0] == b);
                    assert(!b || res[1] == offsetIdx);
                }
            }
        }
        //sw.stop();
        //writeln("Detection: ", sw.peek.usecs / (testIndex.length^^2) / 2);
    }
}
