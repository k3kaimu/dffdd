import std.stdio;
import std.stdio;
import std.getopt;
import std.range;
import std.algorithm;
import std.array;
import std.format;
import std.string;
import std.numeric;
import std.math;
import std.complex;
import std.datetime;

import carbon.stream;

import dffdd.filter.lms;
import dffdd.filter.polynomial;

import dffdd.utils.fft;

enum size_t P = 8;
enum size_t N = 8;
enum size_t Total = 1024*64;
enum real beta = 1E-3;
enum real sampFreq = 400e3;
enum size_t blockSize = 1024;
enum real PLLB = 20;
enum real PLL_K = PLLB / 0.53;
enum real PLLAW = 1.414 * PLL_K;
enum real PLLW2 = PLL_K ^^ 2;
enum real FLLB = 250;
enum real FLLW = FLLB / 0.25;

version = OutputSpectrum;

private F atan(F)(F y, F x)
{
    return x.signbit ? atan2(-y, -x) : atan2(y, x);
}

void main(string[] args)
{
    string sendFilename, recvFilename, outFilename;
    bool bSpeedMode = false;

    getopt(args,
        "sendData", &sendFilename,
        "recvData", &recvFilename,
        "outData", &outFilename,
        "speedMode", &bSpeedMode
    );

    File sendFile = File(sendFilename),
         recvFile = File(recvFilename),
         outFile = File(outFilename, "w");

    recvFile.seek(919 * 8);

    // |x|^^(q), {p:q=p*2+1} = {0: 1, 1: 3, 2: 5, 3: 7}
    auto updater = new PolynomialLMS!((x, p) => x.abs() ^^ (p*2+1), 4, cfloat)
                                     (beta, 1024, 1, 0.5);

    // |x|^^(q-1)*x, {p:q=p*2+1} = {0: 1, 1: 3, 2: 5, 3: 7}
    auto filter = updater.polynomialFilter!((x, p) => x.abs()^^(p*2) * x, 4, cfloat)(N);

    cfloat[] sendBuf = new cfloat[blockSize],
              recvBuf = new cfloat[blockSize],
              outputBuf = new cfloat[blockSize];

    double[] fftResultRecv = new double[blockSize],
             fftResultSIC = new double[blockSize];

    //creal oldPrompt = 0+0i;
    //real oldCarrErr = 0;
    //real carrNCOFreq = 400;
    //immutable DT = blockSize / sampFreq;

    Fft fftObj = new Fft(blockSize);

    //auto nco = lutNCO!(std.math.expi, 1024*16)(carrNCOFreq, DT, 0);
    auto startTime = Clock.currTime;


    foreach(blockIdx; 0 .. Total)
    {
        auto sendGets = sendFile.rawRead(sendBuf);
        auto recvGets = recvFile.rawRead(recvBuf);
        assert(sendGets.length == sendBuf.length);
        assert(recvGets.length == recvBuf.length);

        filter.apply(sendGets, recvGets, outputBuf);

      version(OutputIteration)
      {
        real sum = 0;
        foreach(i, e; outputBuf){
            sum += e.abs()^^2/blockSize;
            if((blockIdx*blockSize+i) % 10 == 0)
            outFile.writefln("%s,%s,", (blockIdx*blockSize+i)/sampFreq, e.abs()^^2);
        }
        outFile.flush();
        //if(blockIdx % 10 == 0)
        //outFile.writefln("%s,%s,", (blockIdx+1)*blockSize/sampFreq, sum);

        if(blockIdx >= 40) throw new Exception("");
      }

      version(OutputSpectrum)
      {
        // fft
        if(!bSpeedMode){
            auto outputSpec = fftObj.fftWithSwap(outputBuf);
            auto recvSpec = fftObj.fftWithSwap(recvBuf);

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] += recvSpec[i].re^^2 + recvSpec[i].im^^2 / 1000;
                fftResultSIC[i] = outputSpec[i].re^^2 + outputSpec[i].im^^2 / 1000;
            }
        }

        if(!bSpeedMode && blockIdx % 1000 == 0){
            foreach(i; 0 .. blockSize)
            {
                auto before = 10*log10(fftResultRecv[i]),
                     after = 10*log10(fftResultSIC[i]);

                outFile.writefln("%s,%s,%s,%s,%s,%s,", blockIdx, i, i*400e3/blockSize-200.0e3, before, after, before - after);
            }
            outFile.flush();

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] = 0;
                fftResultSIC[i] = 0;
            }
        }
      }

        writefln("%s     %s     %s     %s[k samples/s]", blockIdx, recvBuf[$-1].abs, outputBuf[$-1].abs, (blockIdx+1) * blockSize / ((Clock.currTime - startTime).total!"msecs"() / 1000.0L) / 1000.0L);
    }
}
