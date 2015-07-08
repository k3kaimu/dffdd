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
import dffdd.filter.diagonal;
import dffdd.filter.mempoly;

import dffdd.utils.fft;

enum size_t Total = 1024*3;
enum real sampFreq = 400e3;
enum size_t blockSize = 1024;

version = OutputSpectrum;

private F atan(F)(F y, F x)
{
    return x.signbit ? atan2(-y, -x) : atan2(y, x);
}


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

    File sendFile = File(sendFilename),
         recvFile = File(recvFilename),
         outFile = File(outFilename, "w");

    recvFile.seek(seekOffset * 8);

    auto filter2 = {
        auto state = new MemoryPolynomialState!(cfloat, 16,4,2,2,true)(1);

        auto adapter = new LMSAdapter!(typeof(state))(state, 1E-3, 1024, 0.5);

        return new PolynomialFilter!(typeof(state), typeof(adapter))(state, adapter);
    }();

    //auto filter2 = {
    //    auto state = new DiagonalState!(cfloat,
    //                                    (x, xabs, p) => xabs^^(2*p) * x,
    //                                    (xabs, p) => xabs^^(2*p+1),
    //                                    (11+1)/2, 4)();

    //    // set default power
    //    foreach(ref e; state.power) e = 1;

    //    auto adapter = new LMSAdapter!(typeof(state))(state, 1E-4, 1024, 0.5);

    //    return new PolynomialFilter!(typeof(state), typeof(adapter))(state, adapter);
    //}();

    pragma(msg, filter2.state.state.length * filter2.state.state[0].length);

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

        //filter1.apply(sendGets, recvGets, intermBuf);
        filter2.apply(sendGets, recvGets, outputBuf);

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
            auto recvSpec = fftObj.fftWithSwap(recvBuf);

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
                if(abs(freq) < 25e3){
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
        }
      }

        stdout.flush();
    }
}
