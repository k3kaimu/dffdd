module dffdd.dsp.convolution;

import std.complex;
import std.math;
import std.numeric;


/**

*/
C[] convolution(FftObj, C)(FftObj fftObj, in C[] specA, in C[] specB, C[] dst)
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
C[] convolutionPower(FftObj, C)(FftObj fftObj, in C[] specA, in C[] specB, C[] dst)
{
    dst = convolution(fftObj, specA, specB, dst);
    foreach(ref e; dst)
        e = e.re^^2 + e.im^^2;

    return dst;
}
