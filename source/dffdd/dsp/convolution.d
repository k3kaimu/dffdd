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
C[] convolutionPower(FftObj, C)(FftObj fftObj, in FrequencyDomain!(C[]) specA, in FrequencyDomain!(C[]) specB, C[] dst)
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
Nullable!size_t findSignal(FftObj)(FftObj fftObj,
                                   in FrequencyDomain!(C[]) haystackSpec,
                                   in FrequencyDomain!(C[]) needleSpec,
                                   size_t smpPerSym, real thr)
in{
    assert(haystackSpec.length == needleSpec.length);
}
body{
    Nullable!size_t nullV;

    

    //if(haystack.length < needle.length) return nullV;
    //if(needle.length > haystack.length)
    //    swap(haystack, needle);


}