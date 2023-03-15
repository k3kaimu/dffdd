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

alias C = MirComplex!double;

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

    auto isRayleigh = No.isRayleigh;
    enum asWSV = Yes.asWSV;
    enum ModK = 2;              // 1シンボルあたりのビット数（1: BPSK, 2:QPSKなど）
    enum isSpecOP = No.SpecOP;  // 周波数スペクトルプリコーディング
    // auto isRayleigh = No.isRayleigh;

    immutable resultDir = isRayleigh ? "results_Rayleigh" : "results_AWGN";

    import mir.ndslice : map;
    // writeln(makeInterfereMatrix!C(0.5, 1, 8).sliced.map!"a.re");
    
    // if(0)
    static foreach(usePrePost; [Yes.usePrecoding,/* No.usePrecoding*/])
    static foreach(ALGORITHM; ["EP"/*"SVD", "MMSE", "ZF", "EP", "DAMP", "AMP", "GaBP", "QRM-MLD", "Sphere"*/])
    foreach(K; [16, 64, 256])
    {
        if(ALGORITHM == "Sphere" && K > 16)
            continue;

        if(isRayleigh && ALGORITHM == "QRM-MLD" && K == 256)
            continue;

        if(usePrePost && !(ALGORITHM == "EP" || ALGORITHM == "AMP" || ALGORITHM == "DAMP" || ALGORITHM == "GaBP"))
            continue;

        auto NTAPS_ARR = [K/8];
        if(ALGORITHM == "EP" && K == 256)
            NTAPS_ARR = [/*1, 4, 8, 16, */32];

        foreach(NTAPS; isRayleigh ? NTAPS_ARR : [0])
        foreach(ALPHA; [16/16.0, 12/16.0, 10/16.0, /*9/16.0*/ /*16/16.0, 15/16.0, 14/16.0, 13/16.0, 12/16.0, 11/16.0, 10/16.0, 9/16.0, 8/16.0*/]){
            void task(uint K, uint NTAPS, double ALPHA) {
                immutable uint OSRate = 1;
                // immutable uint NFFT = K * OSRate;                   // FFTサイズ
                // immutable uint N = cast(uint)round(K * ALPHA);      // サブキャリア数
                immutable uint N = K;
                immutable uint NFFT = cast(uint)round(K * 1.0 / ALPHA * OSRate);
                immutable uint M = N * OSRate;                      // 観測数
                writefln("%s, K: %s, N: %s, alpha = %s/%s, nTaps=%s", ALGORITHM, NFFT, N, N, NFFT, NTAPS);

                float[] berResults;
                double[] psd;
                // float[] speedResults;

                // enum usePrePost = Yes.usePrecoding;

                // immutable long simTotalBits = 10_000_000;
                // immutable long nChTrial = isRayleigh ? 1_0000 : 1;
                // immutable long simBitsPerTrial = simTotalBits / nChTrial;

                long nChTrial = isRayleigh ? 100 : 1;
                long simBitsPerTrial = isRayleigh ? /*100_000*/ 1_000 : 1_000_000;

                foreach(ebno; iota(0, 21)) {
                    import std.conv : to;

                    size_t numErr, numTot;

                    foreach(iTrial; 0 .. nChTrial) {
                        SimParams!(ModK, ALGORITHM) params;
                        params.usePrecoding = usePrePost;
                        params.asWSV = asWSV;
                        params.isRayleigh = isRayleigh;
                        params.EbN0dB = ebno;
                        params.totalBits = simBitsPerTrial;
                        params.chNTaps = NTAPS;
                        params.chSeed = cast(uint)iTrial;
                        params.M = M;
                        params.N = N;
                        params.NFFT = NFFT;
                        params.Ncp = NFFT / 8;

                        static if(ALGORITHM == "GaBP")
                            auto res = mainImpl!(ModK, "GaBP", usePrePost, asWSV)(M, N, NFFT, ebno, isRayleigh, cast(uint)iTrial, NTAPS * OSRate, simBitsPerTrial, 50, 0.4);
                        else static if(ALGORITHM == "DAMP")
                        {
                            auto res = mainImpl!(ModK, "DAMP", usePrePost, asWSV)(M, N, NFFT, ebno, isRayleigh, cast(uint)iTrial, NTAPS * OSRate, simBitsPerTrial, 50, 1, 1, 1, 1);
                        }
                        else static if(ALGORITHM == "AMP")
                        {
                            auto res = mainImpl!(ModK, "AMP", usePrePost, asWSV)(M, N, NFFT, ebno, isRayleigh, cast(uint)iTrial, NTAPS * OSRate, simBitsPerTrial, 50);
                        }
                        else static if(ALGORITHM == "EP")
                        {
                            auto res = mainImpl(params, 50);
                        }
                        else static if(ALGORITHM == "QRM-MLD")
                        {
                            auto res =  mainImpl!(ModK, "QRM-MLD")(M, N, NFFT, ebno, isRayleigh, cast(uint)iTrial, NTAPS * OSRate, simBitsPerTrial, 100);
                        }
                        else static if(ALGORITHM == "Sphere" || ALGORITHM == "MMSE" || ALGORITHM == "ZF" || ALGORITHM == "SVD" || ALGORITHM == "TDMLD")
                        {
                            auto res =  mainImpl!(ModK, ALGORITHM)(M, N, NFFT, ebno, isRayleigh, cast(uint)iTrial, NTAPS * OSRate, simBitsPerTrial);
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
                    write("%s/ber_%s_%s_%s_%s_%s_%s.json".format(resultDir, ALGORITHM, M, NFFT, NTAPS, usePrePost ? "wP" : "woP", asWSV ? "WSV" : "USV"), toJSON(jv));
                }
                {
                    auto jv = JSONValue(psd.map!(a => JSONValue(a)).array());
                    write("%s/psd_%s_%s_%s_%s_%s_%s.json".format(resultDir, ALGORITHM, M, NFFT, NTAPS, usePrePost ? "wP" : "woP", asWSV ? "WSV" : "USV"), toJSON(jv));
                }
                {
                    auto jv = JSONValue(psd.map!(a => JSONValue(a)).array());
                    write("%s/psd_%s_%s_%s_%s_%s.json".format(resultDir, ALGORITHM, M, K, NTAPS, usePrePost ? "wP" : "woP"), toJSON(jv));
                }
            }

            taskList.append(&task, K, NTAPS, ALPHA);
            // results["kbps_%s_%s_%s_%s".format(ALGORITHM, M, K, usePrePost ? "wP" : "woP")] = JSONValue(speedResults.map!(a => JSONValue(a)).array());
        }
    }

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


struct SimParams(size_t ModK_, string DETECT_)
{
    enum size_t ModK = ModK_;
    enum string DETECT = DETECT_;
    bool usePrecoding = No.usePrecoding;
    bool asWSV = false;
    bool isRayleigh = false;
    uint chSeed = 0;
    uint chNTaps = 0;
    ptrdiff_t totalBits = 10_000_000;
    double EbN0dB = 10;
    uint M = 256;
    uint N = 256;
    uint NFFT = 256;
    uint Ncp = 32;
}


SimResult mainImpl(Params, Args...)(Params params, Args args)
{
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
    immutable size_t block_bits = N * mod.symInputLength * numBlockSyms;
    immutable double ALPHA = N*1.0 / (NFFT / (M*1.0/N));

    size_t total_bits = 0;
    size_t error_bits = 0;

    // 変調行列の生成
    auto modMat = genSEFDMModMatrix!C(N, 1.0 * N / NFFT);

    // OFDMのサブキャリアを圧縮した手法
    if(asWSV) {
        // 電力割当行列
        size_t numUse = cast(uint)round(N * ALPHA);
        auto pamat = matrix!C(M, M, C(0, 0));
        foreach(i; 0 .. numUse) {
            pamat[i, i] = C(sqrt(N * 1.0 / numUse) * sqrt(M * 1.0 / N), 0);
        }

        modMat[] = pamat;

        auto idftM = genDFTMatrix!C(M).H;
        modMat[] = idftM * modMat;
    }

    // チャネル行列
    auto chMat = matrix!C(M, M);
    chMat[] = identity!C(M);

    if(isRayleigh) {
        // i.i.d.レイリーフェージングの巡回行列
        chMat[] = makeCirculantIIDChannelMatrix!C(M, cast(uint)chNTaps, chSeed);
    }

    // チャネル行列の逆行列
    auto invChMat = matrix!C(M, M);
    {
        import dffdd.math.linalg : inv;
        invChMat[] = inv(chMat);
    }

    // auto dftMatN = genDFTMatrix!C(N);
    // auto dftMatM = genDFTMatrix!C(M);

    auto rwMat = matrix!C(M, M);
    auto recvMat = matrix!C(M, N);
    if(usePrecoding) {
        // modMat[] = modMat * makeHaarMatrix!C(N);
        modMat[] = modMat * makeRandomPermutationMatrix!C(N, 0) * genDFTMatrix!C(N);
        // rwMat = makeHaarMatrix!C(M);
        rwMat[] = identity!C(M);
        recvMat[] = rwMat * chMat * modMat;
    } else {
        rwMat[] = identity!C(M);
        recvMat[] = chMat * modMat;
    }

    static if(DETECT == "GaBP")
    {
        auto detector = new GaBPDetector!(C, typeof(mod), GaBPDetectorMode.toBit)(mod, N, M, recvMat.sliced, SIGMA, args);
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
        auto detector = new EPDetector!(C, typeof(mod))(mod, recvMat, SIGMA, args);
    }
    else static if(DETECT == "QRM-MLD")
    {
        rwMat = matrix!C(N, M);
        rwMat[] = modMat.H * invChMat;

        recvMat = matrix!C(N, N);
        recvMat[] = modMat.H * modMat;
        auto detector = makeQRMMLDDetector!C(mod, recvMat, args);
        // assert(OSRate >= 2);
    }
    else static if(DETECT == "Sphere")
    {
        rwMat = matrix!C(N, M);
        rwMat[] = modMat.H * invChMat;

        recvMat = matrix!C(N, N);
        recvMat[] = modMat.H * modMat;
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
    else static if(DETECT == "SVD")
    {
        rwMat = matrix!C(M, M, C(0));
        rwMat[] = identity!C(M);
        {
            import dffdd.math.exprtemplate;
            import kaleidic.lubeck;

            auto svdResult = svd(modMat.sliced);
            modMat[] = modMat * svdResult.vt.lightScope.matrixed.H;
            rwMat[] = svdResult.u.lightScope.matrixed.H;

            auto diag = matrix!C(N, N, C(0));
            foreach(i; 0 .. N)
                diag[i, i] = 1/svdResult.sigma[i];

            modMat[] = modMat * diag;
        }
        rwMat[] = rwMat * invChMat;
        recvMat = matrix!C(M, N, C(0));
        recvMat[] = rwMat * modMat;
        auto detector = makeMMSEDetector!(C, typeof(mod))(mod, recvMat, SIGMA);
    }
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

        if(sum.length != syms.length / N * M)
            sum = vector!C(syms.length / N * M);

        foreach(i; 0 .. syms.length / N) {
            auto s = syms[i*N .. (i+1)*N].sliced.vectored;
            auto v = sum.sliced()[i*M .. (i+1)*M].vectored;
            v[] = modMat * s;
            v[] = chMat * v;
        }

        // 電力スペクトル密度
        foreach(i; 0 .. numBlockSyms) {
            // CP分
            foreach(j; 0 .. Ncp)
                specAnalyzer.put(sum[(i+1)*M - Ncp + j]);

            foreach(j; 0 .. M)
                specAnalyzer.put(sum[i*M + j]);
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


Matrix!(C, Contiguous) genSEFDMModMatrix(C)(uint N, double alpha)
{
    import std.complex : expi;

    auto mat = matrix!C(N, N);
    foreach(i; 0 .. N)
        foreach(j; 0 .. N)
            mat[i, j] = C(expi(-2*PI/N * i * j * alpha).tupleof) / sqrt(N*1.0);
    
    return mat;
}


Matrix!(C, Contiguous) genDFTMatrix(C)(uint NFFT)
{
    import std.complex : expi;

    auto mat = matrix!C(NFFT, NFFT);
    foreach(i; 0 .. NFFT)
        foreach(j; 0 .. NFFT)
            mat[i, j] = C(expi(-2*PI/NFFT * i * j).tupleof) / sqrt(NFFT*1.0);

    return mat;
}


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
