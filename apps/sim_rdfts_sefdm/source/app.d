module app;

import std.algorithm;
import std.complex;
import std.math;
import std.random;
import std.range;
import std.stdio;
import std.traits;

import dffdd.mod;
import dffdd.utils.fft;
import dffdd.mod.sefdm;

import tuthpc.taskqueue;

alias F = double;
alias C = Complex!F;


extern(C) void openblas_set_num_threads(int num_threads);
extern(C) int openblas_get_num_threads();

shared static this()
{
    openblas_set_num_threads(1);
}


void main()
{
    import std.file : mkdirRecurse;

    auto env = defaultJobEnvironment();
    auto taskList = new MultiTaskList!void();

    void run(string dir, uint nFFT, uint nTone, bool isChNorm, uint nChTap, uint nIter, Flag!"isPerfectChEst" isPerfectChEst, size_t numSymChEst, string chEstMethod)
    {
        import std.format;
        import std.json;
        import std.file;
        JSONValue[] berList;

        string filename;
        if(isPerfectChEst)
            filename = format("%s/nFFT%s_isChNorm%s_nTone%s_nChTap%s_nIter%s_perfect.json", dir, nFFT, isChNorm, nTone, nChTap, nIter);
        else
            filename = format("%s/nFFT%s_isChNorm%s_nTone%s_nChTap%s_nIter%s_EstCh%s_%s.json", dir, nFFT, isChNorm, nTone, nChTap, nIter, numSymChEst, chEstMethod);

        if(std.file.exists(filename))
            return;

        writefln!"[nFFT = %s, isChNorm = %s, nTone = %s, nChTap = %s, nIter = %s, Perfect = %s, numSymChEst = %s, chEstMethod = %s]"(nFFT, isChNorm, nTone, nChTap, nIter, isPerfectChEst, numSymChEst, chEstMethod);
        foreach(vEbN0dB; iota(21)) {
            SimParams!(QPSK!C) params;
            params.sefdm.nFFT = nFFT;
            params.sefdm.nTone = nTone;
            params.sefdm.nData = nFFT;
            params.sefdm.nIter = nIter;
            params.sefdm.rndSeed = 0;
            params.vEbN0dB = vEbN0dB;
            params.nChTaps = nChTap;
            params.isChNorm = isChNorm;
            params.isPerfectChEst = isPerfectChEst;
            params.numSymChEst = numSymChEst;
            params.chEstMethod = chEstMethod;
            params.rndSeed = 0;

            auto simResult = runSimImpl(params);
            writefln!"\t%s: %s"(vEbN0dB, simResult.ber);

            berList ~= JSONValue(simResult.ber);

            if(simResult.ber < 1e-7)
                break;
        }
        writeln();

        std.file.write(filename, JSONValue(berList).toString());
    }

    // AWGN
    foreach(nFFT; [64, 128, 256, 512, 1024, 2048]) {
        mkdirRecurse("AWGN");
        foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5]) {
            taskList.append(&run, "AWGN", nFFT, nTone, true, 1, 20, Yes.isPerfectChEst, 1, null);
        }
    }


    // Rayleigh
    foreach(nFFT; [64, 128, 256, 512, 1024, 2048]) {
        mkdirRecurse("Rayleigh");
        foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5])
        foreach(nChTap; [1, 4, 8, 16, 32, 64]) {
            taskList.append(&run, "Rayleigh", nFFT, nTone, false, nChTap, 20, Yes.isPerfectChEst, 1, null);
        }
    }

    foreach(nFFT; [2048]) {
        mkdirRecurse("Rayleigh_Iter");
        foreach(nIter; iota(2, 41, 2)) {
            taskList.append(&run, "Rayleigh_Iter", nFFT, nFFT / 8 * 5, false, 8, nIter, Yes.isPerfectChEst, 1, null);
        }
    }

    foreach(chEstMethod; ["ZF", "MMSE"]) {
        // AWGN, nonPerfectChEst
        foreach(numSymChEst; [0, 1, 2, 4, 8, 16]) {
            foreach(nFFT; [64, 128, 256, 512, 1024, 2048]) {
                mkdirRecurse("AWGN");
                foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5]) {
                    if(numSymChEst == 0)
                        taskList.append(&run, "AWGN", nFFT, nTone, true, 1, 20, Yes.isPerfectChEst, 1, chEstMethod);
                    else
                        taskList.append(&run, "AWGN", nFFT, nTone, true, 1, 20, No.isPerfectChEst, numSymChEst, chEstMethod);
                }
            }
        }

        // Rayleigh, nonPerfectChEst
        foreach(numSymChEst; [0, 1, 2, 4, 8, 16]) {
            foreach(nFFT; [64, 128, 256, 512, 1024, 2048]) {
                mkdirRecurse("Rayleigh");
                foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5])
                foreach(nChTap; [1, 4, 8, 16, 32, 64]) {
                    if(numSymChEst == 0)
                        taskList.append(&run, "Rayleigh", nFFT, nTone, false, nChTap, 20, Yes.isPerfectChEst, 0, chEstMethod);
                    else
                        taskList.append(&run, "Rayleigh", nFFT, nTone, false, nChTap, 20, No.isPerfectChEst, numSymChEst, chEstMethod);
                }
            }
        }
    }


    // foreach(nFFT; [64, 128, 256, 512, 1024, 2048]) {
    //     foreach(nTone; [/*nFFT, nFFT / 8 * 7, nFFT / 8 * 6,*/ nFFT / 8 * 5])
    //     foreach(nChTap; [1, /*nFFT / 32, nFFT / 32 * 2, nFFT / 32 * 3,*/ nFFT / 32 * 4]) {
    //         taskList.append(&run, nFFT, nTone, false, nChTap, 100);
    //     }
    // }


    tuthpc.taskqueue.run(taskList, env);
}


struct SimResult
{
    double ber;
    long errbits;
    long totalbits;
    double modThroughput, demodThroughput;
}


struct RDFTsSEFDMParams(Mod)
{
    Mod mod;
    uint nFFT;
    uint nTone;
    uint rndSeed = 0;
    uint nIter;
    alias nData = nFFT;
}


struct SimParams(Mod)
{
    RDFTsSEFDMParams!Mod sefdm;
    double vEbN0dB;
    size_t maxErrbits = 1000;
    size_t maxTotalbits = 100_000_000;
    size_t minTotalbits = 1_000_000;
    size_t nChTaps = 8;
    bool isChNorm = false;
    bool isPerfectChEst = true;
    string chEstMethod = "ZF";
    size_t numSymChEst = 1;
    size_t rndSeed;
}


SimResult runSimImpl(Mod)(SimParams!Mod params)
{
    import std.datetime.stopwatch;

    immutable double SIGMA2 = 1.0 / (10^^(params.vEbN0dB/10.0) * params.sefdm.mod.symInputLength) * params.sefdm.nFFT / params.sefdm.nData;

    size_t totalbits = 0;
    size_t errbits = 0;

    size_t rndSeed = 0;

    // auto sefdm = new RDFTsSEFDM!Mod(params.sefdm, SIGMA2);
    StopWatch swMod, swDem;

    while(totalbits < params.minTotalbits || errbits < params.maxErrbits )
    {
        if(totalbits > params.maxTotalbits)
            break;

        ++rndSeed;
        auto chFR = makeChannel(params, rndSeed);

        // 送信ビット系列の生成
        Bit[] txbits;
        txbits = genBits(params.sefdm.nData * params.sefdm.mod.symInputLength, rndSeed, txbits);

        // RDFT-s-SEFDM変調
        swMod.start();
        auto txsignal = modulateSEFDM(params.sefdm, txbits);
        swMod.stop();
        auto rxsignal = txsignal.applyChannel(chFR);
        rxsignal = rxsignal.addAWGN(rndSeed, SIGMA2);

        // OFDMによるチャネル推定
        auto estChFR = estimateChannelByOFDM(params, chFR, SIGMA2);

        swDem.start();
        Bit[] rxbits = demodulateSEFDM(params.sefdm, rxsignal, estChFR, SIGMA2);
        // auto rxbits = sefdm.demodulate(txsignal);
        swDem.stop();
        errbits += txbits.zip(rxbits).map!"a[0] != a[1] ? 1 : 0".sum();
        totalbits += txbits.length;
    }


    SimResult result;
    result.totalbits = totalbits;
    result.errbits = errbits;
    result.ber = errbits * 1.0 / totalbits;
    result.modThroughput = totalbits * 1.0 / swMod.peek.total!"usecs" * 1e6;
    result.demodThroughput = totalbits * 1.0 / swDem.peek.total!"usecs" * 1e6;
    return result;
}


C[] modulateSEFDM(Mod)(in ref RDFTsSEFDMParams!Mod params, Bit[] txbits)
in(txbits.length == params.mod.symInputLength * params.nData)
{
    auto sefdm = new RDFTsSEFDM!(Mod, C)(params.mod, params.nFFT, 0, params.nTone, 1, params.nData, null, params.rndSeed);

    C[] dst;
    return sefdm.modulate(txbits, dst);
    // import std;
    // import dffdd.utils.fft;

    // C[] scmod;
    // scmod = params.mod.modulate(txbits, scmod);

    // auto fftw = makeFFTWObject!(TemplateOf!C)(params.nData);

    // // DFTスプレッド
    // fftw.inputs!F[] = scmod[];
    // fftw.fft!F();
    // scmod[] = fftw.outputs!F[];
    // foreach(ref e; scmod) {
    //     e /= sqrt(scmod.length * 1.0);          // FFTをユニタリにする
    //     e /= sqrt(params.nTone * 1.0 / params.nData);  // 1/sqrt(alpha)
    // }

    // // ランダムシャッフル
    // auto rnd = makeRNG(params.rndSeed, "RND_SHUF_SEFDM");
    // scmod = randomShuffle(scmod, rnd);

    // // 帯域圧縮
    // scmod.length = params.nFFT;
    // scmod[params.nTone .. $] = C(0);

    // fftw = makeFFTWObject!(TemplateOf!C)(params.nFFT);
    // fftw.inputs!F[] = scmod;
    // fftw.ifft!F();
    // scmod[] = fftw.outputs!F[];
    // foreach(ref e; scmod) {
    //     e *= sqrt(scmod.length * 1.0);          // IFFTをユニタリにする
    // }

    // return scmod;
}


Bit[] demodulateSEFDM(Mod)(in ref RDFTsSEFDMParams!Mod params, C[] rxsignal, C[] chFR, double N0)
{
    auto sefdm = new RDFTsSEFDM!(Mod, C)(params.mod, params.nFFT, 0, params.nTone, 1, params.nData, null, params.rndSeed, N0, params.nIter);

    import mir.ndslice : sliced;
    import dffdd.math.vector : vectored;
    sefdm.epDetector.channelSingularValues[] = chFR.sliced.vectored;

    Bit[] dst;
    return sefdm.demodulate(rxsignal, dst);
    // import std;
    // import mir.ndslice : sliced;
    // import dffdd.math;
    // import dffdd.detector.ep;
    // import dffdd.utils.fft;

    // // DFTスプレッド
    // auto W1 = dftMatrix!C(params.nData);

    // // ランダムシャッフル
    // auto rnd = makeRNG(params.rndSeed, "RND_SHUF_SEFDM");
    // size_t[] plist = iota(params.nData).map!"cast(size_t)a".array();
    // plist = randomShuffle(plist, rnd);
    // auto P = PermutationMatrix(plist);

    // // 圧縮&電力割当
    // C[] dlist = new C[params.nTone];
    // dlist[] = C(1/sqrt(params.nTone * 1.0 / params.nData));
    // foreach(i; 0 .. params.nTone)
    //     dlist[i] *= chFR[i];

    // // auto W2 = idftMatrix!C(params.nFFT);
    // // EP復調器
    // auto detector = makeSVDEPDetector!C(params.mod, identity!C(params.nTone), dlist.sliced.vectored, P * W1, N0, params.nIter);

    // // 受信信号を周波数領域へ変換
    // {
    //     auto fftw = makeFFTWObject!(TemplateOf!C)(params.nFFT);
    //     fftw.inputs!F[] = rxsignal[];
    //     fftw.fft!F();
    //     rxsignal[] = fftw.outputs!F[];
    //     foreach(ref e; rxsignal)
    //         e /= sqrt(params.nFFT * 1.0);
    // }

    // C[] detected;
    // detected = detector.detect(rxsignal[0 .. params.nTone], detected);

    // Bit[] demodBits;
    // demodBits = params.mod.demodulate(detected, demodBits);
    // return demodBits;
}


Random makeRNG(T)(size_t seed, auto ref T id)
{
    import std;
    return Random(cast(uint) hashOf(forward!id, seed));
}


Bit[] genBits(size_t nbits, size_t seed, return ref Bit[] dst)
{
    import dffdd.utils.binary;
    import std;

    if(dst.length != nbits)
        dst.length = nbits;

    randomBits(makeRNG(seed, "GEN_BITS")).take(nbits).map!(a => Bit(a)).copy(dst);
    return dst;
}


R addAWGN(R)(R received, size_t seed, double SIGMA)
{
    import std;
    import dffdd.blockdiagram.noise;
    auto rnd = BoxMuller!(Random, C)(makeRNG(seed, "GEN_AWGN"));

    foreach(ref e; received) {
        e += rnd.front * sqrt(SIGMA/2);
        rnd.popFront();
    }

    return received;
}


C[] makeChannel(Mod)(in ref SimParams!Mod params, size_t rndSeed)
{
    import dffdd.blockdiagram.noise;
    import dffdd.utils.distribution;

    auto rnd = makeRNG(params.rndSeed + rndSeed, "CHANNEL");
    auto bm = BoxMuller!(Random, C)(rnd);

    C[] impResp = new C[params.sefdm.nFFT];
    impResp[] = C(0);
    foreach(i; 0 .. params.nChTaps) {
        impResp[i] = bm.front / sqrt(params.nChTaps * 2.0);
        bm.popFront();
    }

    auto fftw = makeFFTWObject!(TemplateOf!C)(params.sefdm.nFFT);
    fftw.inputs!F[] = impResp[];
    fftw.fft!F();

    C[] dstFR = new C[params.sefdm.nFFT];
    dstFR[] = fftw.outputs!F[];

    if(params.isChNorm) {
        F sumP = 0;
        foreach(i; 0 .. params.sefdm.nTone)
            sumP += dstFR[i].sqAbs;

        sumP /= params.sefdm.nTone;
        
        foreach(ref e; dstFR)
            e /= sqrt(sumP);
    }
    
    return dstFR;
}


C[] applyChannel(in C[] sig, in C[] chFR)
in(sig.length == chFR.length)
{
    auto fftw = makeFFTWObject!(TemplateOf!C)(sig.length);
    fftw.inputs!F[] = sig[];
    fftw.fft!F();
    foreach(i; 0 .. sig.length)
        fftw.inputs!F[i] = fftw.outputs!F[i] * chFR[i];
    
    fftw.ifft!F();
    foreach(i; 0 .. sig.length)
        fftw.outputs!F[i] = fftw.outputs!F[i];


    return fftw.outputs!F().dup;
}



C[] estimateChannelByOFDM(Mod)(in ref SimParams!Mod params, in C[] chFR, double SIGMA2)
{
    if(params.isPerfectChEst)
        SIGMA2 = 0;
    
    auto sefdm = new RDFTsSEFDM!(Mod, C)(params.sefdm.mod, params.sefdm.nFFT, 0, params.sefdm.nTone, 1, params.sefdm.nData, null, params.sefdm.rndSeed, SIGMA2, params.sefdm.nIter);
    auto ofdm = sefdm.ofdm;

    // チャネル推定用のOFDM信号の生成
    C[] chEstSCtx = new C[params.sefdm.nTone];
    foreach(i; 0 .. params.sefdm.nTone)
        chEstSCtx[i] = i % 2 == 0 ? C(1) : C(-1);
    
    C[] chEstOFDMtx;
    ofdm.modulate(chEstSCtx, chEstOFDMtx);

    C[] estChFR = new C[params.sefdm.nTone];
    estChFR[] = C(0);
    foreach(i; 0 .. params.numSymChEst) {
        auto sig = chEstOFDMtx.dup;
        auto rxsignal = sig.applyChannel(chFR);
        rxsignal = rxsignal.addAWGN(i, SIGMA2);

        C[] chEstSCrx;
        ofdm.demodulate(rxsignal, chEstSCrx);
        // foreach(k; 0 .. params.sefdm.nTone)
        //     estChFR[k] += (chEstSCrx[k] / chEstSCtx[k]) / params.numSymChEst;

        foreach(k; 0 .. params.sefdm.nTone) {
            if(params.chEstMethod == "MMSE")
                estChFR[k] += chEstSCrx[k] / (chEstSCtx[k].sqAbs + SIGMA2) * chEstSCtx[k].conj;
            else if(params.chEstMethod == "ZF")
                estChFR[k] += chEstSCrx[k] / chEstSCtx[k];
            else {
                import std.exception;
                enforce(0, "Unknow channel estimation method");
            }
        }
    }


    foreach(k; 0 .. params.sefdm.nTone) {
        estChFR[k] /= params.numSymChEst;
    }

    return estChFR;

    // if(params.isPerfectChEst) {
    //     return chFR.dup;
    // }

    // C[] estChFR = new C[params.sefdm.nFFT];
    // estChFR[] = C(0);

    // auto fftw = makeFFTWObject!(TemplateOf!C)(params.sefdm.nFFT);

    // // チャネル推定用のOFDMを生成する
    // C[] ofdm;
    // C[] scs;
    // {
    //     foreach(i; 0 .. params.sefdm.nTone)
    //         fftw.inputs!F[i] = (i % 2 == 0) ? 1 : -1;

    //     scs = fftw.inputs!F.dup;
    //     fftw.ifft!F();
    //     foreach(ref e; fftw.outputs!F[]) {
    //         e *= sqrt(params.sefdm.nFFT * 1.0);          // IFFTをユニタリにする
    //     }

    //     ofdm = fftw.outputs!F.dup;
    // }

    // foreach(i; 0 .. params.numSymChEst) {
    //     auto sig = ofdm.dup;
    //     auto rxsignal = sig.applyChannel(chFR);
    //     rxsignal = rxsignal.addAWGN(i, SIGMA2);

    //     fftw.inputs!F[] = rxsignal[];
    //     fftw.fft!F();
    //     foreach(ref e; fftw.outputs!F[]) {
    //         e /= sqrt(params.sefdm.nFFT * 1.0);          // FFTをユニタリにする
    //     }

    //     foreach(j; 0 .. params.sefdm.nTone) {
    //         if(params.chEstMethod == "MMSE")
    //             estChFR[j] += fftw.outputs!F[j] / (scs[j].sqAbs + SIGMA2) * scs[j].conj;
    //         else if(params.chEstMethod == "ZF")
    //             estChFR[j] += fftw.outputs!F[j] / scs[j];
    //         else {
    //             import std.exception;
    //             enforce(0, "Unknow channel estimation method");
    //         }
    //     }
    // }

    // foreach(i; 0 .. params.sefdm.nFFT)
    //     estChFR[i] /= params.numSymChEst;

    // return estChFR;
}
