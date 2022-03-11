module app;

import core.lifetime : forward;

import std.stdio;
import std.range : iota;
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

import dffdd.mod.primitives;
import dffdd.mod.bpsk;
import dffdd.mod.qpsk;
import dffdd.mod.qam;

import dffdd.detector.gabp;
import dffdd.detector.damp;
import dffdd.detector.ep;
import dffdd.detector.qrm_mld;


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
    JSONValue[string] results;
    // {
    //     import std.file : read;
    //     results = parseJSON(cast(const(char)[])read("resultsBER_QPSK.json")).object;
    // }


    static foreach(ALGORITHM; ["EP", "DAMP", "GaBP", "QRM-MLD", "Sphere"])
    foreach(K; [16, 64,/* 256*/])
    {
        if(ALGORITHM == "Sphere" && K > 16)
            continue;

        foreach(ALPHA; [16/16.0, 15/16.0, 14/16.0, 13/16.0, 12/16.0, 11/16.0, 10/16.0]){
            enum ModK = 2;              // 1シンボルあたりのビット数（1: BPSK, 2:QPSKなど）
            immutable uint OSRate = 1;
            immutable uint NFFT = K * OSRate;                   // FFTサイズ
            immutable uint N = cast(uint)round(K * ALPHA);      // サブキャリア数
            immutable uint M = N * OSRate;                      // 観測数
            writefln("%s, K: %s, N: %s, alpha = %s/%s", ALGORITHM, K, N, N, K);

            float[] berResults;
            float[] speedResults;

            foreach(ebno; iota(0, 21)) {
                import std.conv : to;
                // auto res = mainImpl!ModK(M, N, NFFT, ebno, true, args[1].to!uint, args[2].to!uint, args[3].to!double, args[4].to!uint);
                // auto res = mainImpl!(ModK, "DAMP")(M, N, NFFT, ebno, args[1].to!uint, args[2].to!double, args[3].to!double, args[4].to!double, args[5].to!double);
                // auto res = mainImpl!(ModK, "GaBP")(M, N, NFFT, ebno, 10, 0.4);

                static if(ALGORITHM == "GaBP")
                    auto res = mainImpl!(ModK, "GaBP", Yes.usePrecoding)(M, N, NFFT, ebno, 50, 0.4);
                else static if(ALGORITHM == "DAMP")
                {
                    auto res = mainImpl!(ModK, "DAMP", Yes.usePrecoding)(M, N, NFFT, ebno, 50, 1, 1, 1, 1);
                    // auto res = mainImpl!(ModK, "DAMP")(M, N, NFFT, ebno, 10, 0.21392455894027587, 0.048649601825814, 0.7340469614181485, 0.5646985259084107);
                    // auto res = mainImpl!(ModK, "DAMP")(M, N, NFFT, ebno, 20, 0.15853341118873465, 0.023387205571840916, 0.5915693314590689, 0.6479578055604401);
                }
                else static if(ALGORITHM == "EP")
                {
                    auto res = mainImpl!(ModK, "EP", Yes.usePrecoding)(M, N, NFFT, ebno, 50);
                }
                else static if(ALGORITHM == "QRM-MLD")
                {
                    auto res =  mainImpl!(ModK, "QRM-MLD")(M, N, NFFT, ebno, 100);
                }
                else static if(ALGORITHM == "Sphere")
                {
                    auto res =  mainImpl!(ModK, "Sphere")(M, N, NFFT, ebno);
                }
                else static assert(0);

                berResults ~= res.ber;
                speedResults ~= res.kbps;
                writefln!"%s,%s"(ebno, res.ber);
                if(res.ber < 1E-5) {
                    break;
                }
            }

            import std.algorithm : map;
            import std.array : array;
            results["ber_%s_%s_%s".format(ALGORITHM, M, N)] = JSONValue(berResults.map!(a => JSONValue(a)).array());
            results["kbps_%s_%s_%s".format(ALGORITHM, M, N)] = JSONValue(speedResults.map!(a => JSONValue(a)).array());
        }
    }

    import std.file : write;
    auto jv = JSONValue(results);
    write("resultsBER.json", toJSON(jv));
}




size_t lcm(size_t n, size_t m)
{
    return n / gcd(n, m) * m;
}


Tuple!(double, "ber", double, "kbps") mainImpl(size_t ModK, string DETECT, Flag!"usePrecoding" usePrecoding = No.usePrecoding, Args...)(
    uint M, uint N, uint NFFT, double EbN0dB, Args args)
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

    immutable K = mod.symInputLength;

    immutable size_t block_bits = N * K;

    size_t total_bits = 0;
    size_t error_bits = 0;

    // チャネル行列の生成
    // C[][] h = new C[][](N, M);
    auto chMat = matrix!C(M, N);
    foreach(m; 0 .. M) {
        auto s = genChannel!C(N, NFFT, m);
        foreach(n; 0 .. N) {
            chMat[m, n] = s[n].conj;
        }
    }

    // auto dftMatN = genDFTMatrix!C(N);
    // auto dftMatM = genDFTMatrix!C(M);

    static if(usePrecoding)
    {
        chMat[] = chMat * makeHaarMatrix!C(N);
    }
    // chMat[] = chMat * identity!C(N);

    // auto newChMat = matrix!C(N, M);

    // import kaleidic.lubeck : svd;
    // import mir.ndslice : mirtransposed = transposed, mirmap = map;
    // auto svdResult = svd(chMat.sliced);
    // auto weightMat = svdResult.u.mirtransposed.mirmap!(a => conj(a)).slice().matrixed;
    // newChMat[] = weightMat * chMat;

    static if(DETECT == "GaBP")
    {
        // auto rwMat = identity!C(M); // 受信ウェイト行列
        auto rwMat = makeHaarMatrix!C(M);
        auto recvMat = matrix!C(M, N);
        recvMat[] = rwMat * chMat;
        auto detector = new GaBPDetector!(C, typeof(mod), GaBPDetectorMode.toBit)(mod, N, M, recvMat.sliced, SIGMA, args);
    }
    else static if(DETECT == "DAMP")
    {
        // auto rwMat = identity!C(M);
        auto rwMat = makeHaarMatrix!C(M);
        auto recvMat = matrix!C(M, N);
        recvMat[] = rwMat * chMat;
        auto detector = new DAMPDectector!(C, typeof(mod))(mod, recvMat, args);
    }
    else static if(DETECT == "EP")
    {
        // auto rwMat = identity!C(M);
        auto rwMat = makeHaarMatrix!C(M);
        auto recvMat = matrix!C(M, N);
        recvMat[] = rwMat * chMat;
        auto detector = new EPDectector!(C, typeof(mod))(mod, recvMat, SIGMA, args);
    }
    else static if(DETECT == "QRM-MLD")
    {
        auto rwMat = matrix!C(N, M);
        rwMat[] = chMat.H;

        auto recvMat = matrix!C(N, N);
        recvMat[] = rwMat * chMat;
        auto detector = makeQRMMLDDetector!C(mod, recvMat, args);
        // assert(OSRate >= 2);
    }
    else static if(DETECT == "Sphere")
    {
        auto rwMat = matrix!C(N, M);
        rwMat[] = chMat.H;

        auto recvMat = matrix!C(N, N);
        recvMat[] = rwMat * chMat;
        auto detector = makeSphereDetector!C(mod, recvMat, args);
    }
    else
        static assert(0);

    // writeln(chMat.sliced);

    StopWatch sw;
    foreach(_; 0 .. 10_000_000 / block_bits) {
        // writeln(error_bits);
        immutable baseSeed = _;


        // auto bits = genBPSK(M, baseSeed);
        auto bits = genBits(block_bits, baseSeed);
        auto syms = mod.genSymbols(bits);

        // C[] sum = new C[](syms.length / M * N);
        auto sum = vector!C(syms.length / N * M);
        // C[] wsum = new C[](syms.length);
        foreach(i; 0 .. syms.length / N) {
            // mir.blas.gemv(C(1), chMat, syms[i*M .. (i+1)*M].sliced, C(0), sum[i*N .. (i+1)*N].sliced);
            sum.sliced()[i*M .. (i+1)*M].vectored()[] = chMat * syms[i*N .. (i+1)*N].sliced.vectored;
        }

        sum.sliced.addAWGN(baseSeed, SIGMA);


        auto recvY = vector!C(sum.length / rwMat.length!1 * rwMat.length!0);
        foreach(i; 0 .. sum.length / rwMat.length!1) {
            recvY.sliced()[i * rwMat.length!0 .. (i+1) * rwMat.length!0].vectored()[]
                = rwMat * sum.sliced()[i * rwMat.length!1 .. (i+1) * rwMat.length!1].vectored;
        }

        Bit[] decoded;

        sw.start();
        static if(is(detector.OutputElementType == Bit))
        {
            decoded = detector.detect(recvY.sliced.iterator[0 .. recvY.length], decoded);
        }
        else
        {
            C[] decodedSyms;
            decodedSyms = detector.detect(recvY.sliced.iterator[0 .. recvY.length], decodedSyms);
            // writeln(decodedSyms);
            decoded = mod.demodulate(decodedSyms, decoded);
        }
        sw.stop();

        foreach(i; 0 .. decoded.length) {
            total_bits += 1;
            error_bits += bits[i] != decoded[i] ? 1 : 0;
        }

        if(error_bits > 100) {
            break;
        }
    }
    // sw.stop();
    writefln!"--- %s kbps"(total_bits * 1E3 / sw.peek.total!"usecs");

    return typeof(return)(error_bits * 1.0L / total_bits, total_bits * 1E3 / sw.peek.total!"usecs");
}


Random makeRNG(T)(size_t seed, auto ref T id)
{
    return Random(cast(uint) hashOf(forward!id, seed));
}


Bit[] genBits(size_t nbits, size_t seed)
{
    import dffdd.utils.binary;
    import std.range;
    import std.array;
    import std.algorithm;

    return randomBits(makeRNG(seed, "GEN_BITS")).take(nbits).map!(a => Bit(a)).array();
}


C[] genSymbols(Mod)(Mod mod, in Bit[] bits)
{
    C[] dst = new C[](bits.length / mod.symInputLength);
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


Matrix!(C, Contiguous) genDFTMatrix(C)(uint NFFT)
{
    import std.complex : expi;

    auto mat = matrix!C(NFFT, NFFT);
    foreach(i; 0 .. NFFT)
        foreach(j; 0 .. NFFT)
            mat[i, j] = C(expi(2*PI/NFFT * i * j).tupleof) / sqrt(NFFT*1.0);
    
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
