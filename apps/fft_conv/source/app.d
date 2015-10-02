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
import std.typecons;

import carbon.math;
import carbon.stream;

import dffdd.dsp.convolution;
import dffdd.utils.fft;


void main(string[] args)
{
    string sendFilename = "send.dat", recvFilename = "recv.dat"/*, outputFilename = "output.csv"*/;
    size_t blockSize = 1024;
    size_t totalIteration = 1;
    ptrdiff_t offset = 0;

    getopt(args,
        "sendData", &sendFilename,
        "recvData", &recvFilename,
        //"output", &outputFilename,
        "blockSize", &blockSize,
        "totalIteration", &totalIteration,
        "offset", &offset);

    enforce(blockSize.isPowOf2, "Invalid Argument: blockSize is not a power number of 2.");

    Fft fftObj = new Fft(blockSize);

    cfloat[] readBuf = new cfloat[blockSize];
    File sendFile = File(sendFilename);
    File recvFile = File(recvFilename);
    //File outputFile = File(outputFilename, "w");
    Complex!float[] sendData = new Complex!float[blockSize],
                    sendSpec = sendData.dup,
                    recvData = sendData.dup,
                    recvSpec = sendData.dup,
                    rsltData = sendData.dup;

    sendFile.seek(offset*8);

    sendFile.readRawComplex(readBuf, sendData);
    fftObj.fft(sendData, sendSpec);

    foreach(iterIdx; 0 .. totalIteration)
    {
        if(iterIdx == 0)    recvFile.readRawComplex(readBuf, recvData);
        else                recvFile.readRawComplexHalf(readBuf, recvData);
        fftObj.fft(recvData, recvSpec);

        auto res = fftObj.findConvolutionPeak(sendSpec.frequencyDomain, recvSpec.frequencyDomain, rsltData, 40, 10, true);
        if(res.isNull)
            writefln("iteration %s: cannot find.", iterIdx);
        else{
            writefln("iteration %s: find offset: %s -> %s.", iterIdx, iterIdx*blockSize/2+res.get, cast(int)(iterIdx*blockSize/2+res.get) - cast(int)offset);
            break;
        }

        //fftObj.convolutionPower(recvSpec.frequencyDomain, sendSpec.frequencyDomain, rsltData);
        //foreach(i, e; rsltData[0 .. $/2])
        //    outputFile.writefln("%s,%s,", iterIdx*blockSize/2+i, e.re);
    }
}


//void readRawComplex(File file, cfloat[] buf, Complex!float[] output)
//{
//    enforce(file.rawRead(buf).length == buf.length);
//    buf.copyToComplexArray(output);
//}


//void readRawComplexHalf(File file, cfloat[] buf, Complex!float[] output)
//{
//    output[0 .. $/2] = output[$/2 .. $];
//    readRawComplex(file, buf[$/2 .. $], output[$/2 .. $]);
//}


//void copyToComplexArray(in cfloat[] input, Complex!float[] output)
//in{
//    assert(input.length <= output.length);
//}
//body{
//    foreach(i, e; input) output[i] = complex!float(e.re, e.im);
//}


Nullable!size_t findConvolutionPeak(FftObj)(
                    FftObj fftObj,
                    in FrequencyDomain!(Complex!float[]) sendSpec,
                    in FrequencyDomain!(Complex!float[]) recvSpec,
                    Complex!float[] convDst,
                    real dBThreshold = 20,
                    size_t cntThreshold = 4,
                    bool onlyHalf = false)
{
    fftObj.convolutionPower(recvSpec, sendSpec, convDst);

    real maxP = -real.infinity;
    size_t maxIdx;
    foreach(i, ref e; convDst[0 .. onlyHalf ? $/2 : $]){
        auto p = e.re^^2 + e.im^^2;
        e = p;

        if(maxP < p){
            maxP = p;
            maxIdx = i;
        }
    }

    real sum = 0;
    size_t sumCnt;
    foreach(i, e; convDst[0 .. onlyHalf ? $/2 : $]){
        if((max(i, maxIdx) - min(i, maxIdx)) > convDst.length/10){
            sum += e.re;
            ++sumCnt;
        }
    }
    sum /= sumCnt;

    auto snr = maxP / sum;
    auto snrdB = 10*log10(snr);

    Nullable!size_t nullV;

    if(snrdB > 10) writeln("SNR: ", snrdB);

    if(snrdB > dBThreshold)
        return Nullable!size_t(maxIdx);
    else
        return nullV;
}
