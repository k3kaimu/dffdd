import dffdd;

import std;
import gaussian_bp;
import mir.ndslice : slice, sliced, mirmap = map, mirtransposed = transposed;
import mir.blas;
import dffdd.fec;
import dffdd.detector.gabp;
import kaleidic.lubeck;

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
    foreach(N_OBS; [64, 60, 56, 52, 48]){
        enum ModK = 2;              // 1シンボルあたりのビット数（1: BPSK, 2:QPSKなど）
        immutable uint M = 64;              // 入力数
        immutable uint OSRate = 4;
        immutable uint NFFT = M*OSRate;     // FFTサイズ
        immutable uint N = NFFT/64*N_OBS;      // 観測数
        // immutable double SIGMA = M/(10^^(ebno/10.0)) / ModK * OSRate;  // 雑音電力
        writefln("N_OBS: %s, alpha = %s/%s", N_OBS, N_OBS, M);

        foreach(ebno; iota(0, 21) /*[10]*/ ) {
            // auto ber = mainImpl!ModK(M, OSRate, NFFT, N, ebno, true, args[1].to!uint, args[2].to!uint, args[3].to!double, args[4].to!uint);
            auto ber = mainImpl!ModK(M, OSRate, NFFT, N, ebno, true, 5, 1, 0, 20);
            writefln!"%s,%s"(ebno, ber);
            if(ber < 1E-6) {
                break;
            }
            // writeln(log10(ber));
        }
        // writeln();
    }
}



size_t lcm(size_t n, size_t m)
{
    return n / gcd(n, m) * m;
}


double mainImpl(size_t ModK)(uint M, uint OSRate, uint NFFT, uint N, double EbN0dB, bool useFEC = true, uint inner_LDPC_iter = 2, uint inner_GaBP_iter = 2, double dummping = 0.75, uint outer_iter = 20)
// void main()
{
    immutable double SIGMA = M/(10^^(EbN0dB/10.0)) / ModK * OSRate;  // 雑音電力

    Fec!float fec;
    uint FEC_N, FEC_K;
    if(useFEC)
    {
        FEC_N = 648;
        FEC_K = 540;
        fec = new BldpcIeee80211!float(FEC_N, FEC_K, inner_LDPC_iter);
    }
    else
    {
        FEC_N = 540;
        FEC_K = 540;
        fec = new NopBlockFec!(float, x => Bit(x > 1 ? 0 : 1))(FEC_N);
    }

    static if(ModK == 1) {
        auto mod = BPSK!C();
    } else static if(ModK == 2) {
        auto mod = QPSK!C();
    } else {
        auto mod = QAM!C(2^^ModK);
    }

    immutable K = mod.symInputLength;

    size_t block_bits = lcm(M * K, FEC_N) / FEC_N * FEC_K;

    size_t total_bits = 0;
    size_t error_bits = 0;

    // チャネル行列の生成
    // C[][] h = new C[][](N, M);
    auto chMat = slice!C(N, M);
    foreach(m; 0 .. M) {
        auto s = genChannel!C(M, N, NFFT, m);
        foreach(n; 0 .. N) {
            chMat[n][m] = s[n];
        }
    }

    auto svdResult = svd(chMat);
    auto weightMat = svdResult.u.mirtransposed.mirmap!(a => conj(a)).slice();
    auto newChMat = weightMat.mtimes(chMat);
    // auto newChMat = chMat;

    auto gabp = new GaBPDetector!(C, typeof(mod), GaBPDetectorMode.toP0P1)(mod, M, N, newChMat, SIGMA, inner_GaBP_iter, dummping, 4);

    StopWatch sw;
    sw.start();
    foreach(_; 0 .. 10_000_000 / block_bits) {
        // writefln!"%s/%s"(_, 10_000_000 / block_bits);
        // writefln!"iter: %s"(_);
        immutable baseSeed = _;


        // auto bits = genBPSK(M, baseSeed);
        auto bits = genBits(block_bits, baseSeed);
        Bit[] encodedBits;
        fec.encode(bits, encodedBits);
        auto syms = mod.genSymbols(encodedBits);

        C[] sum = new C[](syms.length / M * N);
        // C[] wsum = new C[](syms.length);
        foreach(i; 0 .. syms.length / M) {
            mir.blas.gemv(C(1), newChMat, syms[i*M .. (i+1)*M].sliced, C(0), sum[i*N .. (i+1)*N].sliced);
            // wsum.sliced[i*M .. (i+1)*M] = weightMat.mtimes(sum[i*N .. (i+1)*N].sliced);
        }

        sum.addAWGN(baseSeed, SIGMA);

        // auto wsum = weightMat.mtimes(sum.sliced).array;
        
        Bit[] decodedBits;
        float[] p0p1;
        foreach(i; 0 .. outer_iter) {
            foreach(ref e; p0p1)
                e = log(e);     // to LLR

            float[] gabpP0P1;
            if(p0p1.length) {
                // writefln!"Start GaBP: %s, %s, %s"(sum.length, p0p1.length, gabpP0P1.length);
                gabp.detect(sum, p0p1, gabpP0P1);
                // writefln!"End GaBP: %s, %s, %s"(sum.length, p0p1.length, gabpP0P1.length);
            }
            else {
                // writefln!"Start GaBP: %s, %s"(sum.length, gabpP0P1.length);
                gabp.detect(sum, gabpP0P1);
                // writefln!"End GaBP: %s, %s"(sum.length, gabpP0P1.length);
            }

            if(auto ldpc =  cast(BldpcIeee80211!float)fec) {
                // writefln!"Start LDPC: %s, %s, %s"(gabpP0P1.length, decodedBits.length, p0p1.length);
                ldpc.decode(gabpP0P1, decodedBits, p0p1);
                // writefln!"End LDPC: %s, %s, %s"(gabpP0P1.length, decodedBits.length, p0p1.length);

                // bool errFree = true;
                // foreach(j; 0 .. block_bits / FEC_N)
                //     errFree = errFree && ldpc.checkParity(decodedBits[j*FEC_N .. (j+1)*FEC_N]);

                // if(errFree)
                //     break;
            } else {
                fec.decode(gabpP0P1, decodedBits);

                break;
            }
        }
        
        
        // float[] gabpP0P1;
        // gabp.detect(sum, gabpP0P1);
        // Bit[] decodedBits = new Bit[gabpP0P1.length];
        // fec.decode(gabpP0P1, decodedBits);
        // foreach(i; 0 .. decodedBits.length)
        //     decodedBits[i] = gabpP0P1[i] > 1 ? 0 : 1;
        // decodedBits.length = gabpP0P1.length;
        // foreach(i; 0 .. gabpP0P1.length)
        //     decodedBits[i] = gabpP0P1[i] > 1 ? 0 : 1;

        foreach(i; 0 .. decodedBits.length) {
            total_bits += 1;
            error_bits += bits[i] != decodedBits[i] ? 1 : 0;
        }

        if(error_bits > 100) {
            break;
        }
    }
    sw.stop();
    // writefln!"%s kbps"(total_bits * 1E3 / sw.peek.usecs);

    return error_bits * 1.0L / total_bits;
}


Random makeRNG(T)(size_t seed, auto ref T id)
{
    return Random(cast(uint) hashOf(forward!id, seed));
}


Bit[] genBits(size_t nbits, size_t seed)
{
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
    C[] coefs = new C[N];
    foreach(i, ref e; coefs) {
        e = C(std.complex.expi(2*PI/NFFT * k * i).tupleof);
    }
    
    return coefs;
}



C[] addAWGN(C)(C[] received, size_t seed, double SIGMA)
{
    auto rnd = BoxMuller!(Random, C)(makeRNG(seed, "GEN_AWGN"));

    foreach(ref e; received) {
        e += rnd.front * sqrt(SIGMA/2);
        rnd.popFront();
    }

    return received;
}
