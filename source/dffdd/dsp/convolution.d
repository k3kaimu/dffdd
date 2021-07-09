module dffdd.dsp.convolution;

import std.complex;
import std.math;
import std.numeric;
import std.typecons;
import std.algorithm;

import dffdd.utils.fft;


/**

*/
C[] convolution(FftObj, C)(FftObj fftObj, in C[] specA, in C[] specB, C[] dst)
in{
    assert(specA.length == specB.length);
}
do{
    immutable size = specA.length;
    C[] bufbuf = new C[size];

    if(dst.length < size)
        dst.length = specA.length;

    dst = dst[0 .. size];
    dst[] = specA[];

    foreach(i, ref e; dst)
        e *= -specB[i].conj;

    .ifft!(typeof(C.re))(fftObj, dst, bufbuf);
    dst[] = bufbuf[];

    foreach(i, ref e; dst)
        e /= size;

    return dst;
}


/**

*/
C[] convolutionPower(FftObj, C)(FftObj fftObj,
                                in C[] specA,
                                in C[] specB,
                                C[] dst)
{
    dst = convolution(fftObj, specA, specB, dst);
    foreach(ref e; dst)
        e = e.re^^2 + e.im^^2;

    return dst;
}


struct ConvResult(Index)
{
    Nullable!Index index;
    real snr;

    real snrdB() @property { return 10*log10(snr); }
}



/**
Parameters:
    + smpPerSym := 1シンボルあたりのサンプル数, ただし、OFDMなら1
    + thr := 捕捉しきい値[dB]
*/
ConvResult!size_t findConvolutionPeak(FftObj, C)(
                    FftObj fftObj,
                    in C[] sendSpec,
                    in C[] recvSpec,
                    C[] convDst,
                    real dBThreshold = 20,
                    bool onlyHalf = false)
{
    fftObj.convolutionPower(recvSpec, sendSpec, convDst);

    real maxP = -real.infinity;
    size_t maxIdx;
    foreach(i, ref e; convDst[0 .. onlyHalf ? $/2 : $]){
        auto p = e.re^^2 + e.im^^2;         // ここ2乗必要？
        e = p;

        if(maxP < p){
            maxP = p;
            maxIdx = i;
        }
    }

    real sum = 0;
    size_t sumCnt;
    foreach(i, e; convDst[0 .. onlyHalf ? $/2 : $]){
        if((max(i, maxIdx) - min(i, maxIdx)) > convDst.length/10){
            sum += e.re;
            ++sumCnt;
        }
    }
    sum /= sumCnt;

    auto snr = maxP / sum;
    auto snrdB = 10*log10(snr);

    Nullable!size_t nullV;

    if(snrdB > dBThreshold)
        return typeof(return)(Nullable!size_t(maxIdx), snrdB);
    else
        return typeof(return)(nullV, snrdB);
}