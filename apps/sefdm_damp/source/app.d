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
import dffdd.math.exprtemplate : identity;

import dffdd.mod.primitives;
import dffdd.mod.bpsk;
import dffdd.mod.qpsk;
import dffdd.mod.qam;

import dffdd.detector.gabp;
import dffdd.detector.damp;
import dffdd.detector.qrm_mld;

import mir.ndslice : sliced, slice;

// enum uint M = 64;            // 入力数
// enum uint OSRate = 2;
// enum uint NFFT = M*OSRate;  // FFTサイズ
// enum uint N = NFFT/64*64;     // 観測数
// enum EbN0dB = 10;
// enum ModK = 2;              // 1シンボルあたりのビット数（1: BPSK, 2:QPSKなど）
// enum double SIGMA = M/(10^^(EbN0dB/10.0)) / ModK * OSRate;  // 雑音電力

alias C = MirComplex!double;

void main(string[] args)
{
    JSONValue[string] results;
    static foreach(ALGORITHM; ["DAMP", "GaBP", "QRM-MLD", "Sphere" ])
    foreach(M; [16, 64, 256, /*1024, 4096*/])
    {
        if(ALGORITHM == "Sphere" && M > 16)
            continue;

        foreach(ALPHA; [16/16.0, 15/16.0, 14/16.0, 13/16.0, 12/16.0]){
            enum ModK = 2;              // 1シンボルあたりのビット数（1: BPSK, 2:QPSKなど）
            immutable uint OSRate = 1;
            immutable uint NFFT = M*OSRate;     // FFTサイズ
            immutable uint N = cast(uint)round(NFFT * ALPHA);   // 観測数
            writefln("%s, M: %s, N: %s, alpha = %s/%s", ALGORITHM, M, N, N, M);

            float[] berResults;
            float[] speedResults;

            foreach(ebno; iota(0, 21)) {
                import std.conv : to;
                // auto res = mainImpl!ModK(M, OSRate, NFFT, N, ebno, true, args[1].to!uint, args[2].to!uint, args[3].to!double, args[4].to!uint);
                // auto res = mainImpl!(ModK, "DAMP")(M, OSRate, NFFT, N, ebno, args[1].to!uint, args[2].to!double, args[3].to!double, args[4].to!double, args[5].to!double);
                // auto res = mainImpl!(ModK, "GaBP")(M, OSRate, NFFT, N, ebno, 10, 0.4);

                static if(ALGORITHM == "GaBP")
                    auto res = mainImpl!(ModK, "GaBP")(M, OSRate, NFFT, N, ebno, 10, 0.4);
                else static if(ALGORITHM == "DAMP")
                {
                    auto res = mainImpl!(ModK, "DAMP")(M, OSRate, NFFT, N, ebno, 10, 0.21392455894027587, 0.048649601825814, 0.7340469614181485, 0.5646985259084107);
                    // auto res = mainImpl!(ModK, "DAMP")(M, OSRate, NFFT, N, ebno, 20, 0.15853341118873465, 0.023387205571840916, 0.5915693314590689, 0.6479578055604401);
                }
                else static if(ALGORITHM == "QRM-MLD")
                {
                    auto res =  mainImpl!(ModK, "QRM-MLD")(M, OSRate, NFFT, N, ebno, 10);
                }
                else static if(ALGORITHM == "Sphere")
                {
                    auto res =  mainImpl!(ModK, "Sphere")(M, OSRate, NFFT, N, ebno);
                }
                else static assert(0);

                berResults ~= res.ber;
                speedResults ~= res.kbps;
                writefln!"%s,%s"(ebno, res.ber);
                if(res.ber < 1E-6) {
                    break;
                }
            }

            import std.algorithm : map;
            import std.array : array;
            results["ber_%s_%s_%s".format(ALGORITHM, M, N)] = JSONValue(berResults.map!(a => JSONValue(a)).array());
            results["kbps_%s_%s_%s".format(ALGORITHM, M, N)] = JSONValue(speedResults.map!(a => JSONValue(a)).array());
            // writeln();
        }
    }

    // import std.stdio;
    // File("resultsBER.json", "w").lockingTextWriter.toJSON(jv);
    import std.file : write;
    auto jv = JSONValue(results);
    write("resultsBER.json", toJSON(jv));
}




size_t lcm(size_t n, size_t m)
{
    return n / gcd(n, m) * m;
}


Tuple!(double, "ber", double, "kbps") mainImpl(size_t ModK, string DETECT, Args...)(
    uint M, uint OSRate, uint NFFT, uint N, double EbN0dB, Args args)
{
    // static if(DETECT == "QRM-MLD")
    // {
    //     if(OSRate == 1) {
    //         OSRate = 2;
    //         writefln("QSRate = %s for %s", OSRate, DETECT);
    //     }
    // }

    immutable double SIGMA = M/(10^^(EbN0dB/10.0)) / ModK * OSRate;  // 雑音電力

    static if(ModK == 1) {
        auto mod = BPSK!C();
    } else static if(ModK == 2) {
        auto mod = QPSK!C();
    } else {
        auto mod = QAM!C(2^^ModK);
    }

    immutable K = mod.symInputLength;

    immutable size_t block_bits = M * K;

    size_t total_bits = 0;
    size_t error_bits = 0;

    // チャネル行列の生成
    // C[][] h = new C[][](N, M);
    auto chMat = matrix!C(N, M);
    foreach(m; 0 .. M) {
        auto s = genChannel!C(M, N, NFFT, m);
        foreach(n; 0 .. N) {
            chMat[n, m] = s[n];
        }
    }

    // auto newChMat = matrix!C(N, M);

    // import kaleidic.lubeck : svd;
    // import mir.ndslice : mirtransposed = transposed, mirmap = map;
    // auto svdResult = svd(chMat.sliced);
    // auto weightMat = svdResult.u.mirtransposed.mirmap!(a => conj(a)).slice().matrixed;
    // newChMat[] = weightMat * chMat;

    static if(DETECT == "GaBP")
    {
        auto rwMat = identity!C(N); // 受信ウェイト行列
        auto detector = new GaBPDetector!(C, typeof(mod), GaBPDetectorMode.toBit)(mod, M, N, chMat.sliced, SIGMA, args);
    }
    else static if(DETECT == "DAMP")
    {
        auto rwMat = identity!C(N);
        auto detector = new DAMPDectector!(C, typeof(mod))(mod, chMat, args);
    }
    else static if(DETECT == "QRM-MLD")
    {
        auto rwMat = matrix!C(M, N);
        rwMat[] = chMat.H;

        auto recvMat = matrix!C(M, M);
        recvMat[] = rwMat * chMat;
        auto detector = makeQRMMLDDetector!C(mod, recvMat, args);
        // assert(OSRate >= 2);
    }
    else static if(DETECT == "Sphere")
    {
        auto rwMat = matrix!C(M, N);
        rwMat[] = chMat.H;

        auto recvMat = matrix!C(M, M);
        recvMat[] = rwMat * chMat;
        auto detector = makeSphereDetector!C(mod, recvMat, args);
    }
    else
        static assert(0);

    // writeln(chMat.sliced);

    StopWatch sw;
    foreach(_; 0 .. 10_000_000 / block_bits) {
        immutable baseSeed = _;


        // auto bits = genBPSK(M, baseSeed);
        auto bits = genBits(block_bits, baseSeed);
        auto syms = mod.genSymbols(bits);

        // C[] sum = new C[](syms.length / M * N);
        auto sum = vector!C(syms.length / M * N);
        // C[] wsum = new C[](syms.length);
        foreach(i; 0 .. syms.length / M) {
            // mir.blas.gemv(C(1), chMat, syms[i*M .. (i+1)*M].sliced, C(0), sum[i*N .. (i+1)*N].sliced);
            sum.sliced()[i*N .. (i+1)*N].vectored()[] = chMat * syms[i*M .. (i+1)*M].sliced.vectored;
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


C[] genChannel(C)(uint M, uint N, uint NFFT, uint k)
{
    import std.complex;

    C[] coefs = new C[N];
    foreach(i, ref e; coefs) {
        e = C(std.complex.expi(2*PI/NFFT * k * i).tupleof);
    }
    
    return coefs;
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
