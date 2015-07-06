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
import std.exception;

import carbon.math;
import carbon.stream;

import dffdd.dsp.convolution;
import dffdd.utils.fft;


void main(string[] args)
{
    string sendFilename = "send.dat", recvFilename = "recv.dat", outputFilename = "output.csv";
    size_t blockSize = 1024;
    size_t totalIteration = 1;

    getopt(args,
        "sendData", &sendFilename,
        "recvData", &recvFilename,
        "output", &outputFilename,
        "blockSize", &blockSize,
        "totalIteration", &totalIteration);

    enforce(blockSize.isPowOf2, "Invalid Argument: blockSize is not a power number of 2.");

    Fft fftObj = new Fft(blockSize);

    cfloat[] readBuf = new cfloat[blockSize];
    File sendFile = File(sendFilename);
    File recvFile = File(recvFilename);
    File outputFile = File(outputFilename, "w");
    Complex!float[] sendData = new Complex!float[blockSize],
                    sendSpec = sendData.dup,
                    recvData = sendData.dup,
                    recvSpec = sendData.dup,
                    rsltData = sendData.dup;

    sendFile.readRawComplex(readBuf, sendData);
    fftObj.fft(sendData, sendSpec);

    foreach(iterIdx; 0 .. totalIteration)
    {
        if(iterIdx == 0)    recvFile.readRawComplex(readBuf, recvData);
        else                recvFile.readRawComplexHalf(readBuf, recvData);
        fftObj.fft(recvData, recvSpec);

        fftObj.convolutionPower(recvSpec.frequencyDomain, sendSpec.frequencyDomain, rsltData);
        foreach(i, e; rsltData[0 .. $/2])
            outputFile.writefln("%s,%s,", iterIdx*blockSize/2+i, e.re);
    }
}


void readRawComplex(File file, cfloat[] buf, Complex!float[] output)
{
    enforce(file.rawRead(buf).length == buf.length);
    buf.copyToComplexArray(output);
}


void readRawComplexHalf(File file, cfloat[] buf, Complex!float[] output)
{
    output[0 .. $/2] = output[$/2 .. $];
    readRawComplex(file, buf[$/2 .. $], output[$/2 .. $]);
}


void copyToComplexArray(in cfloat[] input, Complex!float[] output)
in{
    assert(input.length <= output.length);
}
body{
    foreach(i, e; input) output[i] = complex!float(e.re, e.im);
}
