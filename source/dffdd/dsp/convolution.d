module dffdd.dsp.convolution;

import std.complex;
import std.math;
import std.numeric;
import std.typecons;
import std.algorithm;

import dffdd.utils.fft;


/**

*/
C3[] convolution(FftObj, C1, C2, C3)(FftObj fftObj, in FrequencyDomain!(C1[]) specA, in FrequencyDomain!(C2[]) specB, C3[] dst)
in{
    assert(specA.length == specB.length);
}
body{
    alias F = typeof(C1.re);

    immutable size = specA.length;
    // C3[] bufbuf = new C3[size];

    if(dst.length < size)
        dst.length = specA.length;

    dst = dst[0 .. size];
    dst[] = specA[];

    foreach(i, ref e; dst)
        e *= -specB[i].conj;

    // fftObj.inverseFft(dst, bufbuf);
    {
        fftObj.inputs!F[] = dst[];
        fftObj.ifft!F();
    }
    dst[] = fftObj.outputs!F[];

    foreach(i, ref e; dst)
        e /= size;

    return dst;
}


/**

*/
C3[] convolutionPower(FftObj, C1, C2, C3)(FftObj fftObj,
                                in FrequencyDomain!(C1[]) specA,
                                in FrequencyDomain!(C2[]) specB,
                                C3[] dst)
{
    dst = convolution(fftObj, specA, specB, dst);
    foreach(ref e; dst)
        e = e.re^^2 + e.im^^2;

    return dst;
}


/**
Parameters:
    + smpPerSym := 1シンボルあたりのサンプル数, ただし、OFDMなら1
    + thr := 捕捉しきい値[dB]
*/
Nullable!size_t findConvolutionPeak(FftObj, C1, C2, C3)(
                    FftObj fftObj,
                    in FrequencyDomain!(C1[]) sendSpec,
                    in FrequencyDomain!(C2[]) recvSpec,
                    C3[] convDst,
                    real dBThreshold = 20,
                    bool onlyHalf = false)
{
    fftObj.convolutionPower(recvSpec, sendSpec, convDst);

    real maxP = -real.infinity;
    size_t maxIdx;
    foreach(i, ref e; convDst[0 .. onlyHalf ? $/2 : $]){
        auto p = e.re;         // ここ2乗必要？
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
        return Nullable!size_t(maxIdx);
    else
        return nullV;
}