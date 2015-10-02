module dffdd.dsp.convolution;

import std.complex;
import std.math;
import std.numeric;

import dffdd.utils.fft;


/**

*/
C[] convolution(FftObj, C)(FftObj fftObj, in FrequencyDomain!(C[]) specA, in FrequencyDomain!(C[]) specB, C[] dst)
in{
    assert(specA.length == specB.length);
}
body{
    immutable size = specA.length;
    C[] bufbuf = new C[size];

    if(dst.length < size)
        dst.length = specA.length;

    dst = dst[0 .. size];
    dst[] = specA[];

    foreach(i, ref e; dst)
        e *= -specB[i].conj;

    fftObj.inverseFft(dst, bufbuf);
    dst[] = bufbuf[];

    foreach(i, ref e; dst)
        e /= size;

    return dst;
}


/**

*/
C[] convolutionPower(FftObj, C)(FftObj fftObj,
                                in FrequencyDomain!(C[]) specA,
                                in FrequencyDomain!(C[]) specB,
                                C[] dst)
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
Nullable!size_t findConvolutionPeak(FftObj)(
                    FftObj fftObj,
                    in FrequencyDomain!(Complex!float[]) sendSpec,
                    in FrequencyDomain!(Complex!float[]) recvSpec,
                    Complex!float[] convDst,
                    real dBThreshold = 20,
                    bool onlyHalf = false)
{
    fftObj.convolutionPower(recvSpec, sendSpec, convDst);

    real maxP = -real.infinity;
    size_t maxIdx;
    foreach(i, ref e; convDst[0 .. onlyHalf ? $/2 : $]){
        auto p = e.re^^2 + e.im^^2;
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