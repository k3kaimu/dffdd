import std.algorithm;
import std.array;
import std.complex;
import std.datetime;
import std.format;
import std.getopt;
import std.math;
import std.numeric;
import std.path;
import std.range;
import std.stdio;
import std.stdio;
import std.string;

import carbon.stream;

import dffdd.filter.diagonal;
import dffdd.filter.lms;
import dffdd.filter.ls;
import dffdd.filter.mempoly;
import dffdd.filter.polynomial;
import dffdd.utils.fft;

enum size_t Total = 1024*30;
enum real sampFreq = 400e3;
enum size_t blockSize = 1024;

enum size_t N = 4;
enum size_t P = 2;
enum size_t Mb = 0;
enum size_t Mc = 0;
enum bool withDCBias = true;
enum bool withIQImbalance = false;

version = OutputSpectrum;


void main(string[] args)
{
    string sendFilename, recvFilename, outFilename;
    bool bSpeedMode = false, noOutput = false;
    size_t seekOffset;

    getopt(args,
        "sendData", &sendFilename,
        "recvData", &recvFilename,
        "outData", &outFilename,
        "speedMode", &bSpeedMode,
        "noOutput", &noOutput,
        "offset", &seekOffset,
    );

    if(outFilename == null && !noOutput){
        auto currTime = Clock.currTime;
        outFilename = format("output_%s_%d_%d_%d_%d_%s_%s_%s.csv", sendFilename.baseName.stripExtension, N, P, Mb, Mc, withDCBias ? "t" : "f", withIQImbalance ? "t" : "f", currTime.toISOString());
    }

    File sendFile = File(sendFilename),
         recvFile = File(recvFilename),
         outFile/* = File(outFilename, "w")*/;

    if(!noOutput)
        outFile = File(outFilename, "w");

    recvFile.seek(seekOffset * 8);

    auto filter1 = {
        auto state = new MemoryPolynomialState!(cfloat, N, P, Mb, Mc, withDCBias, withIQImbalance)(1);

        //auto adapter = new LSAdapter!(typeof(state))(500);
        auto adapter = new LMSAdapter!(typeof(state))(state, 1E-3, 1024, 0.5);

        return new PolynomialFilter!(typeof(state), typeof(adapter))(state, adapter);
    }();

    pragma(msg, filter1.state.state.length * filter1.state.state[0].length);

    cfloat[] sendBuf = new cfloat[blockSize],
             recvBuf = new cfloat[blockSize],
             intermBuf = new cfloat[blockSize],
             outputBuf = new cfloat[blockSize];

    double[] fftResultRecv = new double[blockSize],
             fftResultSIC = new double[blockSize];

    Fft fftObj = new Fft(blockSize);

    auto startTime = Clock.currTime;


    foreach(blockIdx; 0 .. Total)
    {
        auto sendGets = sendFile.rawRead(sendBuf);
        auto recvGets = recvFile.rawRead(recvBuf);
        assert(sendGets.length == sendBuf.length);
        assert(recvGets.length == recvBuf.length);

        filter1.apply(sendGets, recvGets, outputBuf);

      version(OutputIteration)
      {
        real sum = 0;
        foreach(i, e; outputBuf){
            sum += e.abs()^^2/blockSize;
            if((blockIdx*blockSize+i) % 10 == 0)
            outFile.writefln("%s,%s,", (blockIdx*blockSize+i)/sampFreq, e.abs()^^2);
        }
        outFile.flush();

        if(blockIdx >= 40) throw new Exception("");
      }

      version(OutputSpectrum)
      {
        // fft
        if(!bSpeedMode){
            auto outputSpec = fftObj.fftWithSwap(outputBuf);
            auto recvSpec = fftObj.fftWithSwap(recvGets);

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] += (recvSpec[i].re^^2 + recvSpec[i].im^^2) / 1000;
                fftResultSIC[i] += (outputSpec[i].re^^2 + outputSpec[i].im^^2) / 1000;
            }
        }

        if(!bSpeedMode && blockIdx % 100 == 0){
            real sum = 0; size_t fcnt;
            foreach(i; 0 .. blockSize)
            {
                auto before = 10*log10(fftResultRecv[i]),
                     after = 10*log10(fftResultSIC[i]);

                real freq = i*sampFreq/blockSize-(sampFreq/2);
                if(abs(freq) < 12.5e3){
                    sum += 10.0L ^^(-(before - after)/10);
                    ++fcnt;
                }

                if(!noOutput)
                    outFile.writefln("%s,%s,%s,%s,%s,%s,", blockIdx, i, freq, before, after, before - after);
            }

            writefln("%s,%s,[dB],%s,[k samples/s],", blockIdx, 10*log10(sum / fcnt), (blockIdx+1) * blockSize / ((Clock.currTime - startTime).total!"msecs"() / 1000.0L) / 1000.0L);

            outFile.flush();

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] = 0;
                fftResultSIC[i] = 0;
            }
        }else if(blockIdx % 1000 == 0)
            writefln("%s,%s,[k samples/s],", blockIdx, (blockIdx+1) * blockSize / ((Clock.currTime - startTime).total!"msecs"() / 1000.0L) / 1000.0L);
      }

        stdout.flush();
    }
}
