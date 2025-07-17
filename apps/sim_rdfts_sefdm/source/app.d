module app;

import std.algorithm;
import std.complex;
import std.conv : text;
import std.math;
import std.random;
import std.range;
import std.stdio;
import std.traits;

import dffdd.mod;
import dffdd.utils.fft;
import dffdd.mod.sefdm;
import dffdd.dsp.statistics : makeFilterBankSpectrumAnalyzer;

import tuthpc.taskqueue;

alias F = float;
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

    static void run(alias simImpl, ParamType)(string filename, ParamType[] params, Flag!"shortcut" shortcut = No.shortcut)
    {
        import std.format;
        import std.json;
        import std.file;
        JSONValue[] berList;
        JSONValue[] PAPRs;
        JSONValue[] psd;

        if(std.file.exists(filename ~ "_BER.json"))
            return;

        foreach(p; params) {
            auto simResult = simImpl(p);
            writefln!"\t%s: %s, MOD: %s Mbps, DEMOD: %s Mbps"(p.vEbN0dB, simResult.ber, simResult.modThroughput / 1e6, simResult.demodThroughput / 1e6);

            berList ~= JSONValue(simResult.ber);
            if(p.outputPAPR)
                PAPRs ~= simResult.PAPRs.map!(a => JSONValue(a)).array;

            if(p.outputPSD)
                psd = simResult.psd.map!(a => JSONValue(a)).array;

            if(shortcut && simResult.ber < 1e-7)
                break;
        }
        writeln();

        std.file.write(filename ~ "_BER.json", JSONValue(berList).toString());

        if(PAPRs.length > 0)
            std.file.write(filename ~ "_PAPR.json", JSONValue(PAPRs).toString());

        if(psd.length > 0)
            std.file.write(filename ~ "_PSD.json", JSONValue(psd).toString());
    }

    // // AWGN
    // foreach(nFFT; [64, 128, 256, 512, 1024, 2048]) {
    //     mkdirRecurse("AWGN");
    //     foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5]) {
    //         taskList.append(&run, "AWGN", nFFT, nTone, true, 1, 20, Yes.isPerfectChEst, 1, "Perfect");
    //     }
    // }


    // // Rayleigh
    // foreach(nFFT; [64, 128, 256, 512, 1024, 2048]) {
    //     mkdirRecurse("Rayleigh");
    //     foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5])
    //     foreach(nChTap; [1, 4, 8, 16, 32, 64]) {
    //         taskList.append(&run, "Rayleigh", nFFT, nTone, false, nChTap, 20, Yes.isPerfectChEst, 1, "Perfect");
    //     }
    // }

    // foreach(nFFT; [2048]) {
    //     mkdirRecurse("Rayleigh_Iter");
    //     foreach(nIter; iota(2, 41, 2)) {
    //         taskList.append(&run, "Rayleigh_Iter", nFFT, nFFT / 8 * 5, false, 8, nIter, Yes.isPerfectChEst, 1, "Perfect");
    //     }
    // }

    foreach(chEstMethod; [ChannelEstimationMethod.perfect, ChannelEstimationMethod.ZF]) {
        auto numSymChEstList = (chEstMethod == ChannelEstimationMethod.perfect) ? [1] : [1, 2, 4, 8, 16];

        // AWGN
        if(chEstMethod == ChannelEstimationMethod.perfect)
        foreach(nSpreadIn; [64, 128, 256, 512, 1024, 2048, 4096]) {
            immutable numSymChEst = 1;

            foreach(nFFT; [8, 64, 256, 1024, 4096]) {
                if(nSpreadIn < nFFT)
                    continue;

                immutable string dirName = "AWGN";
                mkdirRecurse(dirName);

                double[] alphaList = [1.0, 7 / 8.0, 6 / 8.0, 5 / 8.0];
                immutable initAlphaListSize = alphaList.length;
                if(nFFT == nSpreadIn) foreach(i; 0 .. min(nFFT, 512) / 2) {
                    alphaList ~= 1.0 - i * 1.0 / min(nFFT, 512);
                }

                foreach(ia, alpha; alphaList) {
                    alias M = QPSK!C;

                    immutable nTone = cast(int)round(nFFT * alpha);
                    immutable nSpreadOut = cast(int)round(nSpreadIn * alpha);
                    SimParams!M[] params;
                    foreach(vEbN0dB; iota(21)) {
                        SimParams!M p;
                        p.sefdm.nFFT = nFFT;
                        p.sefdm.nTone = nTone;
                        p.sefdm.nSpreadIn = nSpreadIn;
                        p.sefdm.nSpreadOut = nSpreadOut;
                        // p.sefdm.nData = nFFT;
                        p.sefdm.nIter = 20;
                        p.sefdm.rndSeed = 0;
                        p.vEbN0dB = vEbN0dB;
                        p.channelType = "AWGN";
                        p.channel = makeAWGNChannel!C(sigma2: 10.0^^(-vEbN0dB / 10.0) / M.symInputLength, seed: 0);
                        p.numSymChEst = numSymChEst;
                        p.chEstMethod = chEstMethod;
                        p.rndSeed = 0;
                        p.update = delegate(ref SimParams!M p, uint seed) {
                            p.rndSeed = seed;
                            p.sefdm.rndSeed = seed;
                            p.channel = makeAWGNChannel!C(sigma2: 10.0^^(-p.vEbN0dB / 10.0) / M.symInputLength, seed: seed);
                        };
                        p.outputPSD = (ia < initAlphaListSize) ? true : false;

                        params ~= p;
                    }

                    string filename = i"$(dirName)/nSpreadIn$(nSpreadIn)_nFFT$(nFFT)_nTone$(nTone)_nIter20_EstCh$(numSymChEst)_$(chEstMethod)".text();
                    taskList.append(&run!(runSimImpl!(C, M), SimParams!M), filename, params, Yes.shortcut);
                }
            }
        }


        // AWGN, change nIter, xaxis = alpha
        if(chEstMethod == ChannelEstimationMethod.perfect)
        foreach(nIter; [1, 5, 10, 20, 40, 80]) {

            immutable string dirName = "AWGN_Change_nIter";
            mkdirRecurse(dirName);

            immutable nSpreadIn = 4096;
            immutable nFFT = nSpreadIn;
            immutable numSymChEst = 1;

            double[] alphaList;
            foreach(i; 0 .. min(nFFT, 512) / 2) {
                alphaList ~= 1.0 - i * 1.0 / min(nFFT, 512);
            }

            foreach(alpha; alphaList) {
                alias M = QPSK!C;

                immutable nTone = cast(int)round(nFFT * alpha);
                immutable nSpreadOut = cast(int)round(nSpreadIn * alpha);
                SimParams!M[] params;
                foreach(vEbN0dB; iota(21)) {
                    SimParams!M p;
                    p.sefdm.nFFT = nFFT;
                    p.sefdm.nTone = nTone;
                    p.sefdm.nSpreadIn = nSpreadIn;
                    p.sefdm.nSpreadOut = nSpreadOut;
                    // p.sefdm.nData = nFFT;
                    p.sefdm.nIter = nIter;
                    p.sefdm.rndSeed = 0;
                    p.vEbN0dB = vEbN0dB;
                    p.channelType = "AWGN";
                    p.channel = makeAWGNChannel!C(sigma2: 10.0^^(-vEbN0dB / 10.0) / M.symInputLength, seed: 0);
                    p.numSymChEst = numSymChEst;
                    p.chEstMethod = chEstMethod;
                    p.rndSeed = 0;
                    p.update = delegate(ref SimParams!M p, uint seed) {
                        p.rndSeed = seed;
                        p.sefdm.rndSeed = seed;
                        p.channel = makeAWGNChannel!C(sigma2: 10.0^^(-p.vEbN0dB / 10.0) / M.symInputLength, seed: seed);
                    };

                    params ~= p;
                }

                string filename = i"$(dirName)/nSpreadIn$(nSpreadIn)_nFFT$(nFFT)_nTone$(nTone)_nIter$(nIter)_EstCh$(numSymChEst)_$(chEstMethod)".text();
                taskList.append(&run!(runSimImpl!(C, M), SimParams!M), filename, params, Yes.shortcut);
            }
        }

        /+

        // AWGN, change shuffle
        if(chEstMethod == ChannelEstimationMethod.perfect)
        {
            immutable nSpreadIn = 4096;
            immutable nFFT = nSpreadIn;

            string dirName = "AWGN_ChangeShuffle";
            mkdirRecurse(dirName);

            immutable double alpha = 5 / 8.0;
            uint[] nShuffleList = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384];

            foreach(nShuffle; nShuffleList) {
                alias M = QPSK!C;

                immutable nTone = cast(int)round(nFFT * alpha);
                immutable nSpreadOut = cast(int)round(nSpreadIn * alpha);
                SimParams!M[] params;
                foreach(vEbN0dB; iota(21)) {
                    SimParams!M p;
                    p.sefdm.nFFT = nFFT;
                    p.sefdm.nTone = nTone;
                    p.sefdm.nSpreadIn = nSpreadIn;
                    p.sefdm.nSpreadOut = nSpreadOut;
                    // p.sefdm.nData = nFFT;
                    p.sefdm.nIter = 20;
                    p.sefdm.rndSeed = 0;
                    p.vEbN0dB = vEbN0dB;
                    p.channelType = "AWGN";
                    p.channel = makeAWGNChannel!C(sigma2: 10.0^^(-vEbN0dB / 10.0) / M.symInputLength, seed: 0);
                    p.numSymChEst = 1;
                    p.chEstMethod = chEstMethod;
                    p.rndSeed = 0;

                    // p.sefdm.shuffleList = iota(nSpreadIn).map!"cast(size_t)a".array;
                    {
                        size_t[] shuffleList = iota(nSpreadIn).map!"cast(size_t)a".array;
                        shuffleNPairs(shuffleList, nShuffle, 0);

                        p.sefdm.shuffleList = shuffleList.dup;
                    }

                    params ~= p;
                }

                string filename = i"$(dirName)/nSpreadIn$(nSpreadIn)_nFFT$(nFFT)_nTone$(nTone)_nIter20_nShuffle$(nShuffle)".text();
                taskList.append(&run!(runSimImpl!(C, M), SimParams!M), filename, params, Yes.shortcut);
            }
        }

        +/

        // AWGN, change block shuffle
        if(chEstMethod == ChannelEstimationMethod.perfect)
        {
            immutable nSpreadIn = 4096;
            immutable nFFT = 4096;

            string dirName = "AWGN_ChangeSwap";
            mkdirRecurse(dirName);

            immutable double alpha = 5 / 8.0;
            uint[] numSwapList = [0, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536];

            foreach(numSwap; numSwapList) {
                alias M = QPSK!C;

                immutable nTone = cast(int)round(nFFT * alpha);
                immutable nSpreadOut = cast(int)round(nSpreadIn * alpha);
                SimParams!M[] params;
                foreach(vEbN0dB; iota(21)) {
                    SimParams!M p;
                    p.sefdm.nFFT = nFFT;
                    p.sefdm.nTone = nTone;
                    p.sefdm.nSpreadIn = nSpreadIn;
                    p.sefdm.nSpreadOut = nSpreadOut;
                    // p.sefdm.nData = nFFT;
                    p.sefdm.nIter = 20;
                    p.sefdm.rndSeed = 0;
                    p.vEbN0dB = vEbN0dB;
                    p.channelType = "AWGN";
                    // p.channel = makeAWGNChannel!C(sigma2: 10.0^^(-vEbN0dB / 10.0) / M.symInputLength, seed: 0);
                    p.numSymChEst = 1;
                    p.chEstMethod = chEstMethod;
                    p.rndSeed = 0;
                    p.update = ((numSwap) => delegate(ref SimParams!M p, uint seed) {
                        p.rndSeed = seed;
                        p.sefdm.rndSeed = seed;
                        p.channel = makeAWGNChannel!C(sigma2: 10.0^^(-p.vEbN0dB / 10.0) / M.symInputLength, seed: seed);

                        size_t[] shuffleList = iota(p.sefdm.nSpreadIn).map!"cast(size_t)a".array;
                        shuffleNPairs(shuffleList, numSwap, seed);
                        p.sefdm.shuffleList = shuffleList.dup;
                    })(numSwap);
                    p.outputPAPR = true;

                    params ~= p;
                }

                string filename = i"$(dirName)/nSpreadIn$(nSpreadIn)_nFFT$(nFFT)_nTone$(nTone)_nIter20_numSwap$(numSwap)".text();
                taskList.append(&run!(runSimImpl!(C, M), SimParams!M), filename, params, Yes.shortcut);
            }
        }


        // Rayleigh
        foreach(numSymChEst; numSymChEstList) {
            foreach(nFFT; [64, 256, 1024, 4096]) {
                mkdirRecurse("Rayleigh");
                foreach(nTone; [nFFT, nFFT / 8 * 7, nFFT / 8 * 6, nFFT / 8 * 5])
                foreach(nChTap; [1, 2, 4, 8, 16, 32, 64, 128])
                {
                    if(nChTap > nFFT / 4)
                        continue;

                    alias M = QPSK!C;
                    immutable nIter = 20;
                    immutable nSpreadIn = 4096;

                    SimParams!M[] params;
                    foreach(vEbN0dB; iota(21)) {
                        SimParams!M p;
                        p.sefdm.nFFT = nFFT;
                        p.sefdm.nTone = nTone;
                        p.sefdm.nSpreadIn = nSpreadIn;
                        p.sefdm.nSpreadOut = nSpreadIn / nFFT * nTone;
                        // p.sefdm.nData = nFFT;
                        p.sefdm.nIter = nIter;
                        p.sefdm.rndSeed = 0;
                        p.vEbN0dB = vEbN0dB;
                        p.channelType = "Rayleigh";
                        p.channel = makeRayleighChannel!C(nTaps: nChTap, sigma2: 10.0^^(-vEbN0dB / 10.0) / M.symInputLength, seed: 0);
                        p.numSymChEst = numSymChEst;
                        p.chEstMethod = chEstMethod;
                        p.rndSeed = 0;
                        p.update = ((nChTap) => delegate(ref SimParams!M p, uint seed) {
                            p.rndSeed = seed;
                            p.sefdm.rndSeed = seed;
                            p.channel = makeRayleighChannel!C(nTaps: nChTap, sigma2: 10.0^^(-p.vEbN0dB / 10.0) / M.symInputLength, seed: seed);
                        })(nChTap);

                        params ~= p;
                    }

                    auto filename = i"Rayleigh/nSpreadIn$(nSpreadIn)_nFFT$(nFFT)_nTone$(nTone)_nChTap$(nChTap)_nIter$(nIter)_EstCh$(numSymChEst)_$(chEstMethod)".text();
                    taskList.append(&run!(runSimImpl!(C, M), SimParams!M), filename, params, Yes.shortcut);
                }
            }
        }


        // Rayleigh, compression factor
        foreach(numSymChEst; numSymChEstList)
        {
            immutable nSpreadIn = 4096;
            immutable nFFT = nSpreadIn;

            double[] alphaList;
            foreach(i; 0 .. min(nFFT, 512) / 2) {
                alphaList ~= 1.0 - i * 1.0 / min(nFFT, 512);
            }

            foreach(alpha; alphaList) {
                alias M = QPSK!C;
                immutable nIter = 20;
                immutable nChTap = 16;
                immutable nTone = cast(int)round(nFFT * alpha);


                SimParams!M[] params;
                foreach(vEbN0dB; iota(21)) {
                    SimParams!M p;
                    p.sefdm.nFFT = nFFT;
                    p.sefdm.nTone = nTone;
                    p.sefdm.nSpreadIn = nSpreadIn;
                    p.sefdm.nSpreadOut = nSpreadIn / nFFT * nTone;
                    // p.sefdm.nData = nFFT;
                    p.sefdm.nIter = nIter;
                    p.sefdm.rndSeed = 0;
                    p.vEbN0dB = vEbN0dB;
                    p.channelType = "Rayleigh";
                    p.channel = makeRayleighChannel!C(nTaps: nChTap, sigma2: 10.0^^(-vEbN0dB / 10.0) / M.symInputLength, seed: 0);
                    p.numSymChEst = numSymChEst;
                    p.chEstMethod = chEstMethod;
                    p.rndSeed = 0;
                    p.update = ((nChTap) => delegate(ref SimParams!M p, uint seed) {
                        p.rndSeed = seed;
                        p.sefdm.rndSeed = seed;
                        p.channel = makeRayleighChannel!C(nTaps: nChTap, sigma2: 10.0^^(-p.vEbN0dB / 10.0) / M.symInputLength, seed: seed);
                    })(nChTap);

                    params ~= p;
                }

                auto filename = i"Rayleigh/nSpreadIn$(nSpreadIn)_nFFT$(nFFT)_nTone$(nTone)_nChTap$(nChTap)_nIter$(nIter)_EstCh$(numSymChEst)_$(chEstMethod)".text();
                taskList.append(&run!(runSimImpl!(C, M), SimParams!M), filename, params, Yes.shortcut);
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
    // foreach(i; 0 .. taskList.length)
    //     taskList[i]();
}



size_t[] makeRandomSelectedIndex(size_t nBlock, size_t nOut, uint rndSeed)
{
    auto rnd = makeRNG(rndSeed, "RND_PUNCTURE");
    size_t[] list = iota(nBlock).map!"cast(size_t)a".array;
    randomShuffle(list, rnd);
    sort(list[0 .. nOut]);
    return list;
}


struct SimResult
{
    double ber;
    long errbits;
    long totalbits;
    double modThroughput, demodThroughput;
    double[] PAPRs;
    float[] psd;
}


struct RDFTsSEFDMParams(Mod)
{
    Mod mod;
    
    uint nFFT;
    uint nTone;
    uint nSpreadIn;
    uint nSpreadOut;
    uint rndSeed = 0;
    uint nIter;
    bool isShuffle = true;
    immutable(size_t)[] shuffleList = null;
    // size_t numShuffle;
    // alias nData = nFFT;
}


enum ChannelEstimationMethod : string
{
    perfect = "Perfect",
    ZF = "ZF",
    MMSE = "MMSE",
}


struct SimParams(Mod)
{
    void delegate(ref SimParams!Mod, uint seed) update;
    RDFTsSEFDMParams!Mod sefdm;
    string channelType;
    WirelessChannel!C channel;
    double vEbN0dB;
    size_t maxErrbits = 1000;
    size_t maxTotalbits = 100_000_000;
    size_t minTotalbits = 1_000_000;
    size_t minTrials = 10000;
    // size_t nChTaps = 8;
    // bool isChNorm = false;
    ChannelEstimationMethod chEstMethod = ChannelEstimationMethod.ZF;
    size_t numSymChEst = 1;
    uint rndSeed;
    bool outputPAPR = false;
    bool outputPSD = false;
}


void shuffleNPairs(size_t[] list, size_t n, uint rndSeed)
{
    auto rnd = makeRNG(rndSeed, "SHUFFLE");
    foreach(_; 0 .. n) {
        auto i = rnd.front % list.length;
        rnd.popFront();
        auto j = rnd.front % list.length;
        rnd.popFront();

        swap(list[i], list[j]);
    }
}


SimResult runSimImpl(C, Mod)(SimParams!Mod params)
{
    import std.datetime.stopwatch;

    // immutable double SIGMA2 = 1.0 / (10^^(params.vEbN0dB/10.0) * params.sefdm.mod.symInputLength) * params.sefdm.nFFT / params.sefdm.nData;

    size_t totalbits = 0;
    size_t errbits = 0;

    uint iTrial = 0;

    // auto sefdm = new RDFTsSEFDM!Mod(params.sefdm, SIGMA2);
    StopWatch swMod, swDem;
    auto specAnalyzer = makeFilterBankSpectrumAnalyzer!C(1_000_000, 256, 256);

    double[] PAPRs;
    while(iTrial < params.minTrials || totalbits < params.minTotalbits || errbits < params.maxErrbits || (params.outputPSD && specAnalyzer.avgCount < 1000))
    {
        if(iTrial > params.minTrials && totalbits > params.maxTotalbits)
            break;

        ++iTrial;

        params.update(params, seed: iTrial);
        // auto chFR = makeChannel(params, iTrial);

        // 送信ビット系列の生成
        Bit[] txbits;
        txbits = genBits(params.sefdm.nSpreadIn * params.sefdm.mod.symInputLength, seed: iTrial, dst: txbits);

        // RDFT-s-SEFDM変調
        swMod.start();
        auto txsignal = modulateSEFDM(params.sefdm, txbits);
        swMod.stop();

        if(params.outputPSD) std.range.put(specAnalyzer, txsignal);
        // auto rxsignal = txsignal.applyChannel(chFR);
        // rxsignal = rxsignal.addAWGN(rndSeed, SIGMA2);
        C[] rxsignal;
        params.channel.apply(txsignal, rxsignal);

        // OFDMによるチャネル推定
        // auto estChFR = estimateChannelByOFDM(params, chFR, SIGMA2);
        auto estChFR = estimateChannelByOFDM!C(params);

        swDem.start();
        Bit[] rxbits = demodulateSEFDM(params.sefdm, rxsignal, estChFR, params.channel.sigma2);
        // auto rxbits = sefdm.demodulate(txsignal);
        swDem.stop();
        errbits += txbits.zip(rxbits).map!"a[0] != a[1] ? 1 : 0".sum();
        totalbits += txbits.length;

        if(params.outputPAPR && PAPRs.length < 1e4) {
            immutable b = params.sefdm.nSpreadIn / params.sefdm.nFFT;
            foreach(i; 0 .. b) {
                auto peak = txsignal[$/b * i .. $/b * (i+1)].map!(a => a.sqAbs).reduce!max;
                auto avg = txsignal[$/b * i .. $/b * (i+1)].map!(a => a.sqAbs).sum / (txsignal.length / b);
                PAPRs ~= peak / avg;
            }
        }
    }


    SimResult result;
    result.totalbits = totalbits;
    result.errbits = errbits;
    result.ber = errbits * 1.0 / totalbits;
    result.modThroughput = totalbits * 1.0 / swMod.peek.total!"usecs" * 1e6;
    result.demodThroughput = totalbits * 1.0 / swDem.peek.total!"usecs" * 1e6;
    if(params.outputPAPR) result.PAPRs = PAPRs.dup;
    if(params.outputPSD) result.psd = specAnalyzer.psd.dup;
    return result;
}


C[] modulateSEFDM(Mod)(in ref RDFTsSEFDMParams!Mod params, Bit[] txbits)
in(txbits.length == params.mod.symInputLength * params.nSpreadIn)
{
    // auto sefdm = new RDFTsSEFDM!(Mod, C)(params.mod, params.nFFT, 0, params.nTone, 1, params.nSpreadIn, params.nSpreadOut, null, params.rndSeed);
    auto sefdm = new RDFTsSEFDM!(Mod, C)(params.mod,
        params.nFFT, nCP: params.nFFT/4, params.nTone, nUpSampling: 1,
        params.nSpreadIn, params.nSpreadOut, shuffle: params.shuffleList, params.rndSeed, N0: 0, params.nIter);

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
    // auto sefdm = new RDFTsSEFDM!(Mod, C)(params.mod, params.nFFT, 0, params.nSpreadIn, params.nSpreadOut, 1, params.nData, null, params.rndSeed, N0, params.nIter);
    auto sefdm = new RDFTsSEFDM!(Mod, C)(params.mod,
        params.nFFT, nCP: params.nFFT/4, params.nTone, nUpSampling: 1,
        params.nSpreadIn, params.nSpreadOut, shuffle: params.shuffleList, params.rndSeed, N0, params.nIter);

    import mir.ndslice : sliced;
    import dffdd.math.vector : vectored;
    // sefdm.epDetector.channelSingularValues[] = chFR.sliced.vectored;
    sefdm.setFrequencyResponse(chFR);

    Bit[] dst;
    return sefdm.demodulate(rxsignal, dst);
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


C[] estimateChannelByOFDM(C, Mod)(ref SimParams!Mod params)
{
    immutable SIGMA2 = params.channel.sigma2;

    auto sefdm = new RDFTsSEFDM!(Mod, C)(params.sefdm.mod,
        params.sefdm.nFFT, nCP: params.sefdm.nFFT/4, params.sefdm.nTone, nUpSampling: 1,
        params.sefdm.nSpreadIn, params.sefdm.nSpreadOut, shuffle: null, params.sefdm.rndSeed, SIGMA2, params.sefdm.nIter);
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

        C[] rxsignal;
        if(params.chEstMethod == ChannelEstimationMethod.perfect) {
            // 雑音なしで完全に推定できると仮定される場合
            params.channel.apply(sig, rxsignal, Yes.ignoreAWGN);
        } else {
            // その他のときは，雑音を含めて受信された信号を使う
            params.channel.apply(sig, rxsignal);
        }

        C[] chEstSCrx;
        ofdm.demodulate(rxsignal, chEstSCrx);

        foreach(k; 0 .. params.sefdm.nTone) {
            if(params.chEstMethod == ChannelEstimationMethod.MMSE)
                estChFR[k] += chEstSCrx[k] / (chEstSCtx[k].sqAbs + SIGMA2) * chEstSCtx[k].conj;
            else if(params.chEstMethod == ChannelEstimationMethod.ZF || params.chEstMethod == ChannelEstimationMethod.perfect)
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
}


class WirelessChannel(C)
{
    import dffdd.filter.state : MultiFIRState;
    import std;
    import dffdd.blockdiagram.noise;


    this(in C[] taps, double sigma2, uint seed)
    {
        _rnd = BoxMuller!(Random, C)(makeRNG(seed, "GEN_AWGN"));
        _scale = sqrt(sigma2 / 2);

        _taps = taps.dup;
        _filt = MultiFIRState!C(1, taps.length);
        foreach(i; 0 .. taps.length)
            _filt.weight[i, 0] = taps[i];
    }


    void apply(in C[] inputs, ref C[] outputs, Flag!"ignoreAWGN" ignoreAWGN = No.ignoreAWGN)
    {
        if(outputs.length != inputs.length)
            outputs.length = inputs.length;

        foreach(i; 0 .. inputs.length) {
            _filt.update(inputs[i]);
            outputs[i] = _filt.output;

            if(!ignoreAWGN) {
                outputs[i] += _rnd.front * _scale;
                _rnd.popFront();
            }
        }
    }


    immutable(C)[] taps() const { return _taps; }
    double sigma2() const { return _scale^^2 * 2; }

  private:
    BoxMuller!(Random, C) _rnd;
    double _scale;
    immutable(C)[] _taps;
    MultiFIRState!C _filt;
}


WirelessChannel!C makeAWGNChannel(C)(double sigma2, uint seed)
{
    C[] taps = [C(1, 0)];
    return new WirelessChannel!C(taps, sigma2, seed);
}


WirelessChannel!C makeRayleighChannel(C)(size_t nTaps, double sigma2, uint seed)
{
    import std;
    import dffdd.blockdiagram.noise;

    auto rnd = makeRNG(seed, "CHANNEL");
    auto bm = BoxMuller!(Random, C)(rnd);

    C[] impResp = new C[nTaps];
    foreach(i; 0 .. nTaps) {
        impResp[i] = bm.front / sqrt(nTaps * 2.0);
        bm.popFront();
    }

    return new WirelessChannel!C(impResp, sigma2, seed);
}
