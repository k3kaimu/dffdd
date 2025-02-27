module app;

import core.lifetime : forward;

import std.stdio;
import std.array;
import std.range;
import std.math;
import std.datetime.stopwatch : StopWatch;
import std.random;
import std.numeric : gcd;
import std.json;
import std.format;
import std.typecons;

import dffdd.math.matrix;
import dffdd.math.vector;
// import dffdd.math.linalg;
import dffdd.math.complex;
import dffdd.math.matrixspecial;
import dffdd.utils.fft;
import dffdd.dsp.statistics : makeFilterBankSpectrumAnalyzer;

import dffdd.mod.primitives;
import dffdd.mod.bpsk;
import dffdd.mod.qpsk;
import dffdd.mod.qam;

import dffdd.detector.gabp;
import dffdd.detector.damp;
import dffdd.detector.amp;
import dffdd.detector.ep;
import dffdd.detector.qrm_mld;
import dffdd.detector.linear;
import dffdd.detector.tridiag_dp;

import tuthpc.taskqueue;


import mir.ndslice : sliced, slice, Contiguous;

// enum uint M = 64;            // 入力数
// enum uint OSRate = 2;
// enum uint NFFT = M*OSRate;  // FFTサイズ
// enum uint N = NFFT/64*64;     // 観測数
// enum EbN0dB = 10;
// enum ModK = 2;              // 1シンボルあたりのビット数（1: BPSK, 2:QPSKなど）
// enum double SIGMA = M/(10^^(EbN0dB/10.0)) / ModK * OSRate;  // 雑音電力

extern(C) void openblas_set_num_threads(int num_threads);
extern(C) int openblas_get_num_threads();

shared static this()
{
    openblas_set_num_threads(1);
}

alias C = StdComplex!double;

void main(string[] args)
{
    auto env = defaultJobEnvironment();
    auto taskList = new MultiTaskList!void();

    // JSONValue[string] results;
    // scope(exit) {
    //     import std.file : write;
    //     auto jv = JSONValue(results);
    //     write("resultsBER_Rayleigh.json", toJSON(jv));
    // }
    // {
    //     import std.file : read;
    //     results = parseJSON(cast(const(char)[])read("resultsBER.json")).object;
    // }

    // auto isRayleigh = No.isRayleigh;
    // enum asOFDM = Yes.asOFDM;
    enum ModK = 2;              // 1シンボルあたりのビット数（1: BPSK, 2:QPSKなど）
    enum isSpecOP = No.SpecOP;  // 周波数スペクトルプリコーディング
    // auto isRayleigh = No.isRayleigh;

    import mir.ndslice : map;
    // writeln(makeInterfereMatrix!C(0.5, 1, 8).sliced.map!"a.re");
    
    // if(0)
    foreach(isRayleigh; [Yes.isRayleigh, No.isRayleigh])
    // foreach(asOFDM; [Yes.asOFDM, No.asOFDM])
    // static foreach(usePrePost; [Yes.usePrecoding, No.usePrecoding])
    static foreach(ALGORITHM; ["EP",/* "SVD", "MMSE", "ZF", "AMP", "GaBP", "QRM-MLD", "Sphere"*/])
    foreach(nFFT; [/*16, 64,*/ /*256,*/ 2048, /*1024, 4096*/])     // FFTサイズ = データ数
    {
        immutable resultDir = isRayleigh ? "results_Rayleigh" : "results_AWGN";

        if(ALGORITHM == "Sphere" && nFFT > 16)
            continue;

        if(isRayleigh && ALGORITHM == "QRM-MLD")
            continue;

        // if(usePrePost && !(ALGORITHM == "EP" || ALGORITHM == "AMP" || ALGORITHM == "DAMP" || ALGORITHM == "GaBP"))
        //     continue;

        // auto NTAPS_ARR = [K/8];
        // if(ALGORITHM == "EP" && K == 256)
        auto NTAPS_ARR = [1, /*4, 8, 16, 32*/];

        foreach(NTAPS; isRayleigh ? NTAPS_ARR : [0])
        foreach(ALPHA; [16/16.0, /*14/16.0, 12/16.0, 10/16.0,*/ /*9/16.0*/ /*16/16.0, 15/16.0, 14/16.0, 13/16.0, 12/16.0, 11/16.0, 10/16.0, 9/16.0, 8/16.0*/]){
            static void task(bool isRayleigh, string resultDir, uint nFFT, uint NTAPS, double ALPHA) {
                immutable uint OSRate = 1;
                // immutable uint NFFT = K * OSRate;                   // FFTサイズ
                // immutable uint N = cast(uint)round(K * ALPHA);      // サブキャリア数
                // immutable uint N = K;
                // immutable uint NFFT = cast(uint)round(K * 1.0 / ALPHA * OSRate);
                // immutable uint M = N * OSRate;                      // 観測数
                // immutable uint N = cast(uint)round(K * ALPHA);
                // immutable uint NFFT = K * OSRate;
                // immutable uint M = N * OSRate;                      // 観測数
                immutable uint nSC = cast(uint)round(nFFT * ALPHA);     // 有効サブキャリア数
                // writefln("%s, K: %s, N: %s, alpha = %s/%s, nTaps=%s", ALGORITHM, K, N, N, NFFT, NTAPS);
                writefln("%s, FFT: %s, Active: %s, alpha = %s/%s", ALGORITHM, nFFT, nSC, nSC, nFFT);

                float[] berResults;
                double[] psd;
                // float[] speedResults;

                // enum usePrePost = Yes.usePrecoding;

                // immutable long simTotalBits = 10_000_000;
                // immutable long nChTrial = isRayleigh ? 1_0000 : 1;
                // immutable long simBitsPerTrial = simTotalBits / nChTrial;

                long nChTrial = isRayleigh ? 100 : 1;
                long simBitsPerTrial = isRayleigh ? /*100_000*/ 1_000 : 1_000_000;

                // bool usePrePost = false;            // ランダムDFTプリコーディングの有無
                // bool asOFDM = false;                // OFDMのサブキャリアに載せるシンボルを圧縮
                static if(ALGORITHM == "EP" || ALGORITHM == "AMP" || ALGORITHM == "GaBP" || ALGORITHM == "DAMP")
                {
                    // usePrePost = true;
                    // asOFDM = true;
                    enum string SEFDMType = "RDFT-s-SEFDM";
                }
                else    // for "SVD", "MMSE", "ZF", "QRM-MLD", "Sphere"
                {
                    // usePrePost = false;
                    // asOFDM = false;
                    enum string SEFDMType = "BasicSEFDM";
                }

                immutable filenameBER = "%s/ber_%s_%s_%s_%s_%s.json".format(resultDir, ALGORITHM, nSC, nFFT, NTAPS, SEFDMType);
                immutable filenamePSD = "%s/psd_%s_%s_%s_%s_%s.json".format(resultDir, ALGORITHM, nSC, nFFT, NTAPS, SEFDMType);

                import std.file : exists;
                if(filenameBER.exists())
                    return;

                foreach(ebno; iota(0, 21)) {
                    import std.conv : to;

                    size_t numErr, numTot;

                    foreach(iTrial; 0 .. nChTrial) {
                        import core.memory : GC;
                        GC.collect();
                        GC.minimize();

                        SimParams!(ModK, SEFDMType, ALGORITHM) params;
                        // params.usePrecoding = usePrePost;
                        // params.asOFDM = asOFDM;
                        params.isRayleigh = isRayleigh;
                        params.EbN0dB = ebno;
                        params.totalBits = simBitsPerTrial;
                        params.chNTaps = NTAPS;
                        params.chSeed = cast(uint)iTrial;
                        params.nSC = nSC;
                        params.nFFT = nFFT;
                        params.nCP = nFFT / 8;

                        static if(ALGORITHM == "GaBP")
                        {
                            // auto res = mainImpl!(ModK, "GaBP", usePrePost, asOFDM)(M, N, NFFT, ebno, isRayleigh, cast(uint)iTrial, NTAPS * OSRate, simBitsPerTrial, 50, 0.4);
                            auto res = mainImpl(params, 50, 0.4);
                        }
                        // else static if(ALGORITHM == "DAMP")
                        // {
                        //     auto res = mainImpl!(ModK, "DAMP", usePrePost, asOFDM)(M, N, NFFT, ebno, isRayleigh, cast(uint)iTrial, NTAPS * OSRate, simBitsPerTrial, 50, 1, 1, 1, 1);
                        // }
                        else static if(ALGORITHM == "AMP")
                        {
                            // auto res = mainImpl!(ModK, "AMP", usePrePost, asOFDM)(M, N, NFFT, ebno, isRayleigh, cast(uint)iTrial, NTAPS * OSRate, simBitsPerTrial, 50);
                            auto res = mainImpl(params, 50);
                        }
                        else static if(ALGORITHM == "EP")
                        {
                            auto res = mainImpl(params, 50);
                        }
                        else static if(ALGORITHM == "QRM-MLD")
                        {
                            // auto res =  mainImpl!(ModK, "QRM-MLD")(M, N, NFFT, ebno, isRayleigh, cast(uint)iTrial, NTAPS * OSRate, simBitsPerTrial, 100);
                            auto res = mainImpl(params, 100);
                        }
                        else static if(ALGORITHM == "Sphere" || ALGORITHM == "MMSE" || ALGORITHM == "ZF" || ALGORITHM == "SVD" || ALGORITHM == "TDMLD")
                        {
                            // auto res =  mainImpl!(ModK, ALGORITHM)(M, N, NFFT, ebno, isRayleigh, cast(uint)iTrial, NTAPS * OSRate, simBitsPerTrial);
                            auto res = mainImpl(params);
                        }
                        else static assert(0);

                        numErr += cast(size_t)round(simBitsPerTrial * res.ber);
                        numTot += simBitsPerTrial;

                        if(iTrial == 0)
                            psd = res.psd;
                    }

                    berResults ~= numErr*1.0 / numTot;
                    // speedResults ~= res.kbps;
                    writefln!"%s,%s"(ebno, berResults[$-1]);
                    if(berResults[$-1] < 1E-6) {
                        break;
                    }
                }

                import std.algorithm : map;
                import std.array : array;
                // results["ber_%s_%s_%s_%s".format(ALGORITHM, M, K, usePrePost ? "wP" : "woP")] = JSONValue(berResults.map!(a => JSONValue(a)).array());
                import std.file : write;
                {
                    auto jv = JSONValue(berResults.map!(a => JSONValue(a)).array());
                    write(filenameBER, toJSON(jv));
                }
                {
                    auto jv = JSONValue(psd.map!(a => JSONValue(a)).array());
                    write(filenamePSD, toJSON(jv));
                }
            }

            if(NTAPS > nFFT / 4)
                continue;

            taskList.append(&task, cast(bool)isRayleigh, resultDir, nFFT, NTAPS, ALPHA);
            // results["kbps_%s_%s_%s_%s".format(ALGORITHM, M, K, usePrePost ? "wP" : "woP")] = JSONValue(speedResults.map!(a => JSONValue(a)).array());
        }
    }

    // writefln!"num of Tasks: %s"(taskList.length);

    // run(taskList, env);
    foreach(i; 0 .. taskList.length)
        taskList[i]();
}




size_t lcm(size_t n, size_t m)
{
    return n / gcd(n, m) * m;
}


struct SimResult
{
    double ber;
    double kbps;
    long errbits;
    long totalbits;
    double[] psd;
}


struct SimParams(size_t ModK_, string SEFDMType_, string DETECT_)
{
    enum size_t ModK = ModK_;
    enum string DETECT = DETECT_;
    enum string SEFDMType = SEFDMType_;
    // bool usePrecoding = No.usePrecoding;
    // bool asOFDM = false;
    bool isRayleigh = false;
    uint chSeed = 0;
    uint chNTaps = 0;
    ptrdiff_t totalBits = 10_000_000;
    double EbN0dB = 10;
    uint nFFT = 256;
    uint nSC = 256;
    uint nCP = 32;
}


SimResult mainImpl(Params, Args...)(Params params, Args args)
{
    import dffdd.math.matrixspecial;

    with(params)
    {

    // static if(DETECT == "QRM-MLD")
    // {
    //     if(OSRate == 1) {
    //         OSRate = 2;
    //         writefln("QSRate = %s for %s", OSRate, DETECT);
    //     }
    // }

    // immutable double SIGMA = 1.0 * M / N /(10^^(EbN0dB/10.0)) / ModK * OSRate;  
    immutable double SIGMA = 1.0 / (10^^(EbN0dB/10.0) * ModK); // 雑音電力

    static if(ModK == 1) {
        auto mod = BPSK!C();
    } else static if(ModK == 2) {
        auto mod = QPSK!C();
    } else {
        auto mod = QAM!C(2^^ModK);
    }

    immutable size_t numBlockSyms = 4;
    immutable size_t block_bits = nFFT * mod.symInputLength * numBlockSyms;
    immutable double ALPHA = nSC * 1.0 / nFFT;

    size_t total_bits = 0;
    size_t error_bits = 0;

    // // 変調行列の生成
    // auto modMat = matrix!C(nFFT, nFFT);
    // if(!asOFDM) {
    //     modMat[] = genSEFDMModMatrix!C(nFFT, ALPHA);
    // } else {  // OFDMのサブキャリアを圧縮した手法
    //     // 電力割当行列
    //     auto pamat = matrix!C(nFFT, nFFT, C(0, 0));
    //     foreach(i; 0 .. nSC) {
    //         pamat[i, i] = C(sqrt(1.0 / ALPHA), 0);
    //     }
    //     modMat[] = pamat;

    //     // auto idftM = genDFTMatrix!C(nFFT).H;
    //     auto idftM = idftMatrix!C(nFFT);
    //     modMat[] = idftM * modMat;
    // }


    static if(SEFDMType == "RDFT-s-SEFDM")
    {
        // auto modMat = (){
        auto dft_ = dftMatrix!C(nFFT);
        auto rndperm_ = (){
            auto perm = new size_t[nFFT];
            foreach(i; 0 .. nFFT)
                perm[i] = i;
            
            Random rnd;
            rnd.seed(0);
            perm = randomShuffle(perm, rnd);
            return PermutationMatrix(perm);
        }();
        C[] palist;
        foreach(i; 0 .. nFFT)
            palist ~= C( i < nSC ? sqrt(1.0 / ALPHA) : 0 );

        auto pamat_ = diag(palist);
        auto idft_ = idftMatrix!C(nFFT);

        auto modMat = idft_ * pamat_ * rndperm_ * dft_;
        // }();

        enum bool hasIDFT_SVD = true;
        auto modMatSVs = palist;            // 変調行列の特異値のリスト
        auto precodeMat = rndperm_ * dft_;
    }
    else static if(SEFDMType == "BasicSEFDM")
    {
        auto modMat = genSEFDMModMatrix!C(nFFT, ALPHA);
        enum bool hasIDFT_SVD = false;
    }
    else static assert(0, "SEFDMType '" ~ SEFDMType ~ "' is not supported.");

    // // チャネル行列
    // auto chMat = matrix!C(nFFT, nFFT);
    // chMat[] = identity!C(nFFT);

    // if(isRayleigh) {
    //     // i.i.d.レイリーフェージングの巡回行列
    //     chMat[] = makeCirculantIIDChannelMatrix!C(nFFT, cast(uint)chNTaps, chSeed);
    // }


    alias ChannelMatrixType = typeof(idftMatrix!C(nFFT) * diag([C(1)]) * dftMatrix!C(nFFT));
    C[] freqResp;
    ChannelMatrixType chMat;
    ChannelMatrixType invChMat;
    if(isRayleigh) {
        freqResp = makeRayleighFreqResp!C(nFFT, nSC, cast(uint)chNTaps, chSeed);
        chMat = idftMatrix!C(nFFT) * diag(freqResp) * dftMatrix!C(nFFT);

        auto invFreqResp = freqResp.dup;
        foreach(ref e; invFreqResp) e = 1/e;
        invChMat = idftMatrix!C(nFFT) * diag(invFreqResp) * dftMatrix!C(nFFT);
    } else {
        foreach(i; 0 .. nFFT) freqResp ~= C(1);
        
        chMat = idftMatrix!C(nFFT) * diag(freqResp) * dftMatrix!C(nFFT);
        invChMat = chMat;
    }


    // pragma(msg, typeof(chMat * modMat));


    // // チャネル行列の逆行列
    // auto invChMat = matrix!C(nFFT, nFFT);
    // {
    //     import dffdd.math.linalg : inv;
    //     invChMat[] = inv(chMat);
    // }

    // auto dftMatN = genDFTMatrix!C(N);
    // auto dftMatM = genDFTMatrix!C(M);

    // auto rwMat = matrix!C(nFFT, nFFT);
    // auto recvMat = matrix!C(nFFT, nFFT);
    // if(usePrecoding) {
    //     // modMat[] = modMat * makeHaarMatrix!C(N);
    //     modMat[] = modMat * makeRandomPermutationMatrix!C(nFFT, 0) * genDFTMatrix!C(nFFT);
    //     // rwMat = makeHaarMatrix!C(M);
    //     rwMat[] = identity!C(nFFT);
    //     recvMat[] = rwMat * chMat * modMat;
    // } else {
    //     rwMat[] = identity!C(nFFT);
    //     recvMat[] = chMat * modMat;
    // }

    static if(DETECT == "GaBP")
    {
        auto detector = new GaBPDetector!(C, typeof(mod), GaBPDetectorMode.toBit)(mod, nFFT, nFFT, recvMat.sliced, SIGMA, args);
    }
    else static if(DETECT == "DAMP")
    {
        auto detector = new DAMPDectector!(C, typeof(mod))(mod, recvMat, args);
    }
    else static if(DETECT == "AMP")
    {
        auto detector = new AMPDetector!(C, typeof(mod))(mod, recvMat, SIGMA, args);
    }
    else static if(DETECT == "EP")
    {
        // auto recvMat = matrix!C(nFFT, nFFT);
        auto recvMat = chMat * modMat;
        auto rwMat = dftMatrix!C(nFFT);
        // auto rwMat = identity!C(nFFT);
        // auto detector = new EPDetector!(C, typeof(mod))(mod, recvMat, SIGMA, args);
        C[] recvSVs = modMatSVs.dup;
        foreach(i; 0 .. nFFT)
            recvSVs[i] *= freqResp[i];

        auto detector = makeSVDEPDetector!C(mod, identity!C(nFFT), recvSVs.sliced.vectored, precodeMat, SIGMA, args);
    }
    else static if(DETECT == "QRM-MLD")
    {
        // rwMat = matrix!C(nFFT, nFFT);
        // rwMat[] = modMat.H * invChMat;

        // recvMat = matrix!C(nFFT, nFFT);
        auto recvMat = modMat.H * modMat;
        auto detector = makeQRMMLDDetector!C(mod, recvMat, args);
        // assert(OSRate >= 2);
    }
    else static if(DETECT == "Sphere")
    {
        // rwMat = matrix!C(nFFT, nFFT);
        // rwMat[] = modMat.H * invChMat;

        // recvMat = matrix!C(nFFT, nFFT);
        auto recvMat = modMat.H * modMat;
        auto detector = makeSphereDetector!C(mod, recvMat, args);
    }
    else static if(DETECT == "MMSE")
    {
        auto detector = makeMMSEDetector!(C, typeof(mod))(mod, recvMat, SIGMA);
    }
    else static if(DETECT == "ZF")
    {
        auto detector = makeZFDetector!(C, typeof(mod))(mod, recvMat);
    }
    // else static if(DETECT == "SVD")
    // {
    //     rwMat = matrix!C(nFFT, nFFT, C(0));
    //     rwMat[] = identity!C(nFFT);
    //     {
    //         import dffdd.math.exprtemplate;
    //         import kaleidic.lubeck;

    //         auto svdResult = svd(modMat.sliced);
    //         modMat[] = modMat * svdResult.vt.lightScope.matrixed.H;
    //         rwMat[] = svdResult.u.lightScope.matrixed.H;

    //         auto diag = matrix!C(nFFT, nFFT, C(0));
    //         foreach(i; 0 .. nFFT)
    //             diag[i, i] = 1/svdResult.sigma[i];

    //         modMat[] = modMat * diag;
    //     }
    //     rwMat[] = rwMat * invChMat;
    //     recvMat = matrix!C(nFFT, nFFT, C(0));
    //     recvMat[] = rwMat * modMat;
    //     auto detector = makeMMSEDetector!(C, typeof(mod))(mod, recvMat, SIGMA);
    // }
    else 
        static assert(0);

    // writeln(modMat.sliced);

    // size_t sumPSDN = 0;
    // double[] psd = new double[(M + Ncp) * numBlockSyms];
    // psd[] = 0;
    // auto fftw = makeFFTWObject!MirComplex((M + Ncp) * numBlockSyms);
    // scope(exit) {
    //     psd[] /= sumPSDN;
    // }

    // FilterBankSpectrumAnalyze
    auto specAnalyzer = makeFilterBankSpectrumAnalyzer!C(1_000_000, 256, 256);

    Bit[] bits;
    C[] syms;
    Vector!(C, Contiguous) sum;
    Vector!(C, Contiguous) recvY;
    Bit[] decoded;
    C[] decodedSyms;

    StopWatch sw;
    foreach(_; 0 .. (params.totalBits == -1 ? 10_000_000 : params.totalBits) / block_bits + 1) {
        // writeln(error_bits);
        immutable baseSeed = _;

        bits = genBits(block_bits, baseSeed, bits);
        syms = mod.genSymbols(bits, syms);

        if(sum.length != syms.length)
            sum = vector!C(syms.length);

        foreach(i; 0 .. syms.length / nFFT) {
            auto s = syms[i*nFFT .. (i+1)*nFFT].sliced.vectored;
            auto v = sum.sliced()[i*nFFT .. (i+1)*nFFT].vectored;
            v[] = modMat * s;
            v[] = chMat * v;
        }

        // 電力スペクトル密度
        foreach(i; 0 .. numBlockSyms) {
            // CP分
            foreach(j; 0 .. nCP)
                specAnalyzer.put(sum[(i+1)*nFFT - nCP + j]);

            foreach(j; 0 .. nFFT)
                specAnalyzer.put(sum[i*nFFT + j]);
        }

        sum.sliced.addAWGN(baseSeed, SIGMA);

        if(recvY.length != sum.length / rwMat.length!1 * rwMat.length!0)
            recvY = vector!C(sum.length / rwMat.length!1 * rwMat.length!0);

        foreach(i; 0 .. sum.length / rwMat.length!1) {
            recvY.sliced()[i * rwMat.length!0 .. (i+1) * rwMat.length!0].vectored()[]
                = rwMat * sum.sliced()[i * rwMat.length!1 .. (i+1) * rwMat.length!1].vectored;
        }

        // Bit[] decoded;

        sw.start();
        static if(is(detector.OutputElementType == Bit))
        {
            decoded = detector.detect(recvY.sliced.iterator[0 .. recvY.length], decoded);
        }
        else
        {
            decodedSyms = detector.detect(recvY.sliced.iterator[0 .. recvY.length], decodedSyms);
            // writeln(decodedSyms);
            decoded = mod.demodulate(decodedSyms, decoded);
        }
        sw.stop();

        foreach(i; 0 .. decoded.length) {
            total_bits += 1;
            error_bits += bits[i] != decoded[i] ? 1 : 0;
        }

        // writefln("total_bits: %s, error_bits: %s", total_bits, error_bits);

        if(/*simTotalBits == -1 && */error_bits > 100) {
            break;
        }
    }
    // sw.stop();
    // writefln!"--- %s kbps"(total_bits * 1E3 / sw.peek.total!"usecs");

    // return typeof(return)(error_bits * 1.0L / total_bits, total_bits * 1E3 / sw.peek.total!"usecs", error_bits, total_bits, psd);
    SimResult result;
    result.ber = error_bits * 1.0L / total_bits;
    result.errbits = error_bits;
    result.totalbits = total_bits;
    result.kbps = total_bits * 1E3 / sw.peek.total!"usecs";
    result.psd = specAnalyzer.psd.dup;
    // writefln!"alpha=%s, %s (dB): %s"(ALPHA, EbN0dB, result.ber);

    return result;
    }
}


Random makeRNG(T)(size_t seed, auto ref T id)
{
    return Random(cast(uint) hashOf(forward!id, seed));
}


Bit[] genBits(size_t nbits, size_t seed, return ref Bit[] dst)
{
    import dffdd.utils.binary;
    import std.range;
    import std.array;
    import std.algorithm;

    if(dst.length != nbits)
        dst.length = nbits;

    randomBits(makeRNG(seed, "GEN_BITS")).take(nbits).map!(a => Bit(a)).copy(dst);
    return dst;
}


C[] genSymbols(Mod)(Mod mod, in Bit[] bits, return ref C[] dst)
{
    if(dst.length != bits.length / mod.symInputLength)
        dst.length = bits.length / mod.symInputLength;

    mod.modulate(bits, dst);
    return dst;
}


C[] genChannel(C)(uint N, uint NFFT, uint k)
{
    import std.complex : expi;

    C[] coefs = new C[N];
    foreach(i, ref e; coefs) {
        e = C(expi(2*PI/NFFT * k * i).tupleof) / sqrt(N*1.0);
    }
    
    return coefs;
}


// Matrix!(C, Contiguous) genSEFDMModMatrix(C)(uint N, double alpha)
// {
//     import std.complex : expi;

//     auto mat = matrix!C(N, N);
//     foreach(i; 0 .. N)
//         foreach(j; 0 .. N)
//             mat[i, j] = C(expi(-2*PI/N * i * j * alpha).tupleof) / sqrt(N*1.0);
    
//     return mat;
// }

auto genSEFDMModMatrix(C)(uint N, double alpha)
{
    import dffdd.math.matrixspecial;

    immutable uint M = cast(uint)round(N / alpha);

    C[] ones;
    foreach(i; 0 .. N)
        ones ~= C(1);

    auto diag1 = diag(M, N, ones);      // 入力N, 出力Mで，出力に(M-N)のゼロを追加する行列
    auto idft = idftMatrix!C(M, M);     // M x MのIDFT行列
    auto diag2 = diag(N, M, ones);      // 入力M, 出力Nで，出力の(M-N)個の要素を取り除く
    
    return diag2 * idft * diag1;
}


// Matrix!(C, Contiguous) genDFTMatrix(C)(uint NFFT)
// {
//     import std.complex : expi;

//     auto mat = matrix!C(NFFT, NFFT);
//     foreach(i; 0 .. NFFT)
//         foreach(j; 0 .. NFFT)
//             mat[i, j] = C(expi(-2*PI/NFFT * i * j).tupleof) / sqrt(NFFT*1.0);

//     return mat;
// }


Matrix!(C, Contiguous) makeHaarMatrix(C)(uint N)
{
    import dffdd.math.linalg : qrDecomp;

    C[] arr = new C[N * N];
    arr[] = C(0);
    addAWGN(arr, N, 1);

    auto matQ = matrix!C(N, N);
    auto matR = matrix!C(N, N);

    arr.sliced(N, N).matrixed.qrDecomp(matQ, matR);
    return matQ;
}



R addAWGN(R)(R received, size_t seed, double SIGMA)
{
    import dffdd.blockdiagram.noise;
    auto rnd = BoxMuller!(Random, C)(makeRNG(seed, "GEN_AWGN"));

    foreach(ref e; received) {
        e += rnd.front * sqrt(SIGMA/2);
        rnd.popFront();
    }

    return received;
}


Matrix!(C, Contiguous) makeCirculantIIDChannelMatrix(C)(size_t N, size_t nTaps, uint seed)
{
    import std.range : take;
    import std.algorithm : map;
    import std.array : array;
    import dffdd.blockdiagram.noise;
    auto rnd = BoxMuller!(Random, C)(makeRNG(seed, "GEN_CHANNEL"));

    immutable coef = sqrt(1.0/nTaps/2);
    C[] taps = rnd.take(nTaps).map!(a => a * C(coef) ).array();
    C[] row = new C[N];

    row[0] = taps[0];
    row[$ - (nTaps-1) .. $] = taps[1 .. $];

    auto mat = matrix!C(N, N);
    foreach(i; 0 .. N) {
        mat.sliced()[i, 0 .. $] = row;
        {
            auto last = row[N-1];
            foreach_reverse(j; 1 .. N)
                row[j] = row[j-1];

            row[0] = last;
        }
    }

    return mat;
}


// C[] makeRayleighFreqResp(C, alias ComplexTemplate = StdComplex)(size_t nFFT, size_t nTaps, uint seed)
// in(nFFT >= nTaps)
// {
//     import std.range : take;
//     import std.algorithm : map;
//     import std.array : array;
//     import dffdd.blockdiagram.noise;

//     auto rnd = BoxMuller!(Random, C)(makeRNG(seed, "GEN_CHANNEL"));

//     immutable c = sqrt(1.0/nTaps/2);
//     C[] taps = rnd.take(nTaps).map!(a => a * C(c) ).array();

//     foreach(i; 0 .. nFFT - nTaps)
//         taps ~= C(0);

//     alias R = typeof(C.init.re);

//     import dffdd.utils.fft;
//     auto fftw = makeFFTWObject!ComplexTemplate(nFFT);
//     fftw.inputs!R[] = taps;
//     fftw.fft!R();
//     auto coefs = fftw.outputs!R().dup;

//     // R p = 0;
//     // foreach(ref e; coefs) p += e.sqAbs;

//     // immutable normCoef = sqrt(p / nFFT);
//     // foreach(ref e; coefs) {
//     //     e /= normCoef;
//     // }

//     return coefs;
// }


C[] makeRayleighFreqResp(C)(size_t nFFT, size_t nSC, size_t nTaps, uint seed)
{
    import dffdd.blockdiagram.noise;
    import dffdd.utils.distribution;
    import std.traits;
    alias F = typeof(C.init.re);

    auto rnd = makeRNG(seed, "GEN_CHANNEL");
    auto bm = BoxMuller!(Random, C)(rnd);

    C[] impResp = new C[nFFT];
    impResp[] = C(0);
    foreach(i; 0 .. nTaps) {
        impResp[i] = bm.front / sqrt(nTaps * 2.0);
        bm.popFront();
    }

    auto fftw = makeFFTWObject!(TemplateOf!C)(nFFT);
    fftw.inputs!F[] = impResp[];
    fftw.fft!F();

    C[] dstFR = new C[nFFT];
    dstFR[] = fftw.outputs!F[];

    // if(params.isChNorm) {
    //     F sumP = 0;
    //     foreach(i; 0 .. nSC)
    //         sumP += dstFR[i].sqAbs;

    //     sumP /= nSC;
        
    //     foreach(ref e; dstFR)
    //         e /= sqrt(sumP);
    // }
    
    return dstFR;
}



Matrix!(C, Contiguous) makeRandomPermutationMatrix(C)(size_t N, uint seed)
{
    auto rnd = makeRNG(seed, "GEN_Perm");

    auto mat = matrix!C(N, N);
    mat[] = C(0);

    size_t[] idxs = iota(N).array;
    idxs = randomShuffle(idxs, rnd);

    foreach(i; 0 .. N) {
        mat[i, idxs[i]] = 1;
    }

    return mat;
}


Matrix!(C, Contiguous) makePartialSelectMatrix(C)(size_t N, size_t M)
{
    immutable C scale = C(sqrt(N / M * 1.0));
    auto mat = matrix!C(N, N);
    mat[] = C(0);
    foreach(i; 0 .. M)
        mat[i, i] = scale;
    
    return mat;
}
