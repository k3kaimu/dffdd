module dffdd.filterbank.protofilter;

import std.algorithm;
import std.complex;
import std.math;
import std.mathspecial;
import std.array;


import dffdd.window;
import dffdd.math.math;
import dffdd.utils.fft;

import mir.ndslice;


/**
h(n) = sin(wc * (n - 0.5 N)) / {pi * (n - 0.5 N)}
where:
    wc = 2 * pi * cutoff
*/
F[] desingSincPulse(F = real)(size_t ntaps, real cutoff) pure nothrow @safe
{
    F[] sincpulse = new F[ntaps];
    foreach(i; 0 .. ntaps) {
        immutable F x = (i - ntaps/F(2)) * cutoff;
        sincpulse[i] = sinc(x) * cutoff;
    }

    sincpulse[] /= sum(sincpulse);

    return sincpulse;
}

unittest 
{
    // numpyResults = arr / arr.sum()
    // where:
    // arr = np.sinc((np.arange(12) - 6.0)/2)/2
    real[] numpyResults = [
        1.85320633e-17,   6.05303122e-02,  -1.85320633e-17,
        -1.00883854e-01,   1.85320633e-17,   3.02651561e-01,
         4.75403961e-01,   3.02651561e-01,   1.85320633e-17,
        -1.00883854e-01,  -1.85320633e-17,   6.05303122e-02
    ];

    real[] coefs = desingSincPulse!real(12, 0.5);
    assert(std.algorithm.equal!((a,b) => isClose(a-b, 0, 0, 1e-8))(coefs, numpyResults));
}


/**
See also: https://jp.mathworks.com/help/dsp/ref/designmultiratefir.html
*/
F[] designKaiserSincPulse(F)(size_t nchannel, size_t polyphaseTaps, real attdB, real cutoffCoef = 1, real beta = real.nan) pure nothrow @safe
{
    immutable size_t N = nchannel * polyphaseTaps;
    auto sincpulse = desingSincPulse!F(N, cutoffCoef/nchannel);

    if(beta.isNaN)
        beta = kaiserBetaFromStopbandAttdB(attdB);

    F[] kw = kaiser!F(N, beta);
    foreach(i; 0 .. N)
        sincpulse[i] *= kw[i];

    sincpulse[] /= sum(sincpulse);
    
    return sincpulse;
}

unittest
{
    // testResults = coefs / coefs.sum()
    // where:
    // beta = 4.55126
    // coefs = kw * pluse
    // kw = numpy.kaiser(12, beta)
    // pulse =  np.sinc((np.arange(12) - 6.0)/2)/2
    real[] numpyResults = [
        1.08476232e-18,   1.27547671e-02,  -8.07660847e-18,
        -6.99242733e-02,   1.70465256e-17,   3.18657108e-01,
         5.00545414e-01,   2.78390890e-01,   1.28448806e-17,
        -4.39670087e-02,  -3.90502117e-18,   3.54310261e-03
    ];

    numpyResults[] /= sum(numpyResults);

    real[] res = designKaiserSincPulse!real(2, 6, 50, 1);
    assert(std.algorithm.equal!((a, b) => isClose(a-b, 0, 0, 1e-8))(res, numpyResults));
}


/**
0~PIまでの所望の周波数応答を達成するプロトタイプフィルタを設計します
*/
F[] designPrototypeFromFreqResponse(F = real)(in F[] freqResp)
in{
    import carbon.math : isPowOf2;
    assert(freqResp.length.isPowOf2);
}
do {
    import std.algorithm : stdmap = map;

    immutable size_t N = freqResp.length * 2;

    auto fftw = globalBankOf!makeFFTWObject[N];
    auto ips = fftw.inputs!F;
    foreach(i, ref e; ips[0 .. N/2])
        e = freqResp[i];

    foreach(i; 1 .. N/2)
        ips[$-i] = ips[i];

    ips[$/2] = 0;

    fftw.ifft!F();
    auto dst = fftw.outputs!F.stdmap!"a.re".array();
    swapHalf(dst);
    dst[] /= sum(dst);

    return dst;
}



/**
See also: https://jp.mathworks.com/matlabcentral/fileexchange/15813-near-perfect-reconstruction-polyphase-filterbank
*/
F[] designRootRaisedERF(F = real)(size_t nchannel, size_t polyphaseTaps, F k = F.nan, F cutoffCoef = 1)
{
    alias C = Complex!F;

    if(k.isNaN){
        switch(polyphaseTaps){
            case 8:     k = 4.853; break;
            case 10:    k = 4.775; break;
            case 12:    k = 5.257; break;
            case 14:    k = 5.736; break;
            case 16:    k = 5.856; break;
            case 18:    k = 7.037; break;
            case 20:    k = 6.499; break;
            case 22:    k = 6.483; break;
            case 24:    k = 7.410; break;
            case 26:    k = 7.022; break;
            case 28:    k = 7.097; break;
            case 30:    k = 7.755; break;
            case 32:    k = 7.452; break;
            case 48:    k = 8.522; break;
            case 64:    k = 9.396; break;
            case 96:    k = 10.785; break;
            case 128:   k = 11.5; break; 
            case 192:   k = 11.5; break;
            case 256:   k = 11.5; break;
            default:    k = 8;
        }
    }

    return designPrototypeFromFreqDomainWindow!F(nchannel, polyphaseTaps, (F x) => erfc(k*x)*0.5, cutoffCoef);
}

unittest
{
    import std.range;
    import mir.ndslice.dynamic;

    auto coeff = designRootRaisedERF(32, 32, real.nan, 0.5);

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


    import std.algorithm : stdmap = map;
    auto p = sqrt(testResult.stdmap!"a^^2".sum());
    foreach(a, b; lockstep(coeff, testResult.sliced(32, 32).transposed.flattened)){
        assert(abs(a - b)/p < 1E-6);
    }
}


/**
See also: https://en.wikipedia.org/wiki/Root-raised-cosine_filter
*/
F[] designRootRaisedCosine(F = real)(size_t nchannel, size_t polyPhaseTaps, F beta = 0.5, F cutoffCoef = 1)
{
    immutable size_t N = nchannel * polyPhaseTaps;
    immutable F M = F(nchannel) / cutoffCoef;

    
    F[] coefs = new F[N];
    foreach(i, ref e; coefs){
        immutable real n = F(i) - F(N)/2;
        if(n == 0){
            e = 1/M * (1 + beta * (4/PI - 1));
        }else if(abs(n) == M/(4*beta)){
            e = beta / M / SQRT2 * ((1 + 2/PI) * sin(PI/4/beta) + (1 - 2/PI) * cos(PI/4/beta));
        }else{
            immutable F num = sin(PI * n / M * (1 - beta)) + 4 * beta * n / M * cos(PI * n / M * (1 + beta));
            immutable F den =  PI * n / M * (1 - (4*beta*n/M)^^2);
            e = num / den / M;
        }
    }

    coefs[] /= sum(coefs);

    return coefs;
}

unittest
{
    import std.stdio;
    auto coefs = designRootRaisedCosine!real(4, 8, 0.5);

    auto testResults = [
        -0.00252179, -0.000953692, 0.00267476, 0.00410571, 0.000756536, -0.00410571, -0.00374466, 0.00386018,
        0.0105915, 0.00386018, -0.0187233, -0.0391411, -0.0264787, 0.0391411, 0.144401, 0.243191,
        0.283651, 0.243191, 0.144401, 0.0391411, -0.0264787, -0.0391411, -0.0187233, 0.00386018,
        0.0105915, 0.00386018, -0.00374466, -0.00410571, 0.000756536, 0.00410571, 0.00267476, -0.000953692
    ];

    foreach(i; 0 .. coefs.length)
        assert(isClose(coefs[i] - testResults[i], 0, 0, 1e-5));
}



/**
*/
F[] designKaiserBesselDrived(F = real)(size_t nchannel, size_t polyphaseTaps, F beta = 14, F cutoffCoef = 1)
in{
    assert(polyphaseTaps % 2 == 0);
}
do {
    immutable F i0beta = besselI0(beta);

    /**
    int [0 .. x] (1-4t^2)^k dt を計算します
    */
    static
    F int1m4t(F x, F k)
    {
        F dst = x * hypgeom21!F(0.5, -k, 3.0/2.0L, 4.0L*x^^2);
        return dst;
    }

    /*
    w(x) = I0(beta * sqrt(1-4x^2)) / I0(beta) について
    int [0 .. x] w(x) dxを計算します
    */
    static
    F integralKaiser(F beta, F x, F i0beta)
    {
        if(x > 0.5)
            return integralKaiser(beta, 0.5, i0beta);
        else if(x < -0.5)
            return integralKaiser(beta, -0.5, i0beta);

        F sum = 0;
        F r = 1;
        size_t k = 0;
        do {
            ++k;
            r *= (beta/k/2)^^2;
            sum += r * int1m4t(x, k);
        } while(abs(r) > F.epsilon);

        return sum / i0beta;
    }

    // int[-0.5 .. 0.5] w(x) dx
    immutable F intAllKaiser = integralKaiser(beta, 0.5, i0beta) * 2;

    // w(x) = I0(beta * sqrt(1-4x^2)) / I0(beta)を積分して周波数領域上で窓関数を構築
    return
    designPrototypeFromFreqDomainWindow!F(
        nchannel, polyphaseTaps,
        (real x) => 0.5 - integralKaiser(beta, x, i0beta) / intAllKaiser,
        cutoffCoef);
}


/**
wI(x) = int[x -> infty] w(x) dx = int[-\infty -> -x] w(x) dx という窓関数を積分した関数wI(x)からナイキストフィルタを構築します．
ただし，wI(infty)=1でありw(x)は偶関数でw(0)で最大値を取るような窓関数を前提にアルゴリズムが構築されています．
また，一つのチャネル（バンド）に対してw(x)はx=-0.5 ~ 0.5が相当し，全バンドではx=-nchannel/2 ~ nchannel/2まで動きます．
*/
F[] designPrototypeFromFreqDomainWindow(F = real, Fn)(size_t nchannel, size_t polyphaseTaps, scope Fn wI, F cutoffCoef = 1.0)
{
    immutable size_t N = nchannel * polyphaseTaps;
    F[] freqResp = new F[N / 2];
    foreach(i; 0 .. N/2) {
        immutable F freq = F(i) / N / cutoffCoef;
        freqResp[i] = sqrt(wI(freq*nchannel - 0.5) - wI(freq*nchannel + 0.5));
    }

    return designPrototypeFromFreqResponse(freqResp);
}
