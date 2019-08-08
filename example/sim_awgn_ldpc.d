/++ dub.json:
{
    "name": "sim_awgn_ldpc",
    "dependencies": {
        "dffdd": { "path": ".." }
    }
}
+/
import dffdd.fec.fec;
import dffdd.fec.ldpc;
import dffdd.mod.bpsk;
import dffdd.mod.qpsk;
import dffdd.mod.qam;
import dffdd.mod.primitives;

import std.algorithm;
import std.array;
import std.complex;
import std.datetime;
import std.math;
import std.parallelism;
import std.random;
import std.range;
import std.stdio;



void main()
{
    alias C = Complex!float;

    immutable
        N = 648,
        K = 324,
        M = N - K,
        R = K * 1.0 / N;

    auto ldpc = new BldpcIeee80211!float(N, K);
    auto mod = BPSK!(Complex!float)();

    immutable Nc = mod.symInputLength;

    immutable
        MaxTotalBlocks = 10_000,        // 1万ブロックまで伝送する
        MaxErrorBlocks  = 100;          // 100ブロック誤った時点で終わり

    auto ebn0_dBList = iota(0, 5.5, 0.25).array;
    double[] berList = new double[ebn0_dBList.length];
    double[] blerList = new double[ebn0_dBList.length];

    foreach(i_ebn0, ebn0_dB; ebn0_dBList) {
        immutable double ebn0 = 10.0^^(ebn0_dB / 10);
        immutable double snr = ebn0 * K / N * Nc;
        immutable double N0 = 1 / snr;
        immutable double noiseAmp = sqrt(N0);

        auto p0p1calc = p0p1Calculator!float(mod);
        auto decoder = ldpc.makeDecoder();

        Bit[] info = new Bit[K];
        C[] noise = new C[N / Nc];
        C[] received = new C[N / Nc];

        size_t totalBlocks = 0,
               errorBlocks = 0,
               totalBits = 0,
               errorBits = 0;

        Random rnd;
        rnd.seed(0);

        StopWatch sw;
        while(totalBlocks < MaxTotalBlocks && errorBlocks < MaxErrorBlocks) {
            enum size_t P = 8;
            float[] p0p1, llr;
            Bit[] coded;
            Bit[] decoded = new Bit[K*P];
            
            foreach(i; 0 .. P) {
                rnd.makeBits(info);

                Bit[] coded_part;
                ldpc.encode(info, coded_part);
                coded ~= coded_part;

                C[] sym;
                mod.modulate(coded_part, sym);
                rnd.makeNoise(noise);

                foreach(j; 0 .. sym.length)
                    sym[j] += noise[j] * noiseAmp;

                float[] p0p1_part;
                p0p1calc.computeP0P1(sym, p0p1_part, N0);
                p0p1 ~= p0p1_part;
            }

            sw.start();
            decoder.decode(p0p1, decoded);
            sw.stop();

            foreach(i; 0 .. P) {
                if(coded[i*N .. i*N + K] != decoded[i*K .. (i+1)*K])
                    ++errorBlocks;

                foreach(j; 0 .. K)
                    if(coded[i*N + j] != decoded[i*K + j])
                        ++errorBits;
            }

            totalBits += K * P;
            totalBlocks += P;
        }

        berList[i_ebn0] = errorBits * 1.0 / totalBits;
        blerList[i_ebn0] =  errorBlocks * 1.0 / totalBlocks;
        writefln!"EbN0 = %s [dB], DecodeThroughput = %s [kbps]"(
            ebn0_dB,
            totalBits / (sw.peek.usecs / 1.0E6) / 1E3
        );
    }

    foreach(i, e; ebn0_dBList) {
        writefln!"EbN0 = %s [dB], BER = %s, BLER = %s"(e, berList[i], blerList[i]);
    }
}


void makeNoise(C)(ref Random rnd, C[] dst)
{
    foreach(i; 0 .. dst.length) {
        immutable 
            x = uniform01(rnd),
            y = uniform01(rnd);

        immutable r = sqrt(-log(x));
        dst[i] = r * std.complex.expi(2*PI*y);
    }
}


void makeBits(ref Random rnd, Bit[] dst)
{
    foreach(ref e; dst) {
        e = rnd.front % 2;
        rnd.popFront();
    }
}
