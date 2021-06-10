import std.algorithm;
import std.array;
import std.complex;
import std.exception;
import std.file;
import std.format;
import std.getopt;
import std.json;
import std.math;
import std.numeric;
import std.range;
import std.stdio;
import std.string;
import std.typecons;
import std.path;

import carbon.math;
import carbon.stream;

import dffdd.dsp.convolution;
import dffdd.utils.fft;
import dffdd.utils.jsvisitor;


void main(string[] args)
{
    static struct Visitor
    {
        Visitor onDirectory(VisitorData data)
        {
            return this;
        }


        void onFile(VisitorData data)
        {
            if(!data.txFileName.exists)
                return;

            auto res = implMain(data.txFileName, data.rxFileName);
            if(!res.index.isNull){
                writeln("found: ", buildPath(data.parentKeys), ", idx=", res.index.get, ", snr[dB]=", res.snrdB);
                data.entry["offset"] = res.index.get;
            }
            else{
                writeln("not found: ", buildPath(data.parentKeys));
            }
        }
    }

    Visitor v;

    JSONValue jv = args[1].readText().parseJSON();
    visitJSON(jv, v);

    std.file.write("offset_result.json", jv.toPrettyString());
}


ConvResult!ptrdiff_t implMain(string txDataFileName, string rxDataFileName)
{
    ptrdiff_t offsetTX = 200_000;
    size_t blockSize = 8192;
    size_t totalIteration = 100;

    return peakSearch(txDataFileName, rxDataFileName, blockSize, totalIteration, offsetTX);
}


ConvResult!ptrdiff_t peakSearch(string txDataFileName, string rxDataFileName, size_t blockSize, size_t totalIteration, ptrdiff_t offsetTX)
{
    enforce(blockSize.isPowOf2, "Invalid Argument: blockSize is not a power number of 2.");

    Fft fftObj = new Fft(blockSize);

    cfloat[] readBuf = new cfloat[blockSize];
    File txFile = File(txDataFileName);
    File rxFile = File(rxDataFileName);
    //File outputFile = File(outputFilename, "w");
    Complex!float[] sendData = new Complex!float[blockSize],
                    sendSpec = sendData.dup,
                    recvData = sendData.dup,
                    recvSpec = sendData.dup,
                    rsltData = sendData.dup;

    txFile.seek((1_000_000 + offsetTX) * 8);
    rxFile.seek(1_000_000 * 8);

    txFile.readRawComplex(readBuf, sendData);
    fftObj.fft(sendData, sendSpec);

    foreach(iterIdx; 0 .. totalIteration)
    {
        if(iterIdx == 0)    rxFile.readRawComplex(readBuf, recvData);
        else                rxFile.readRawComplexHalf(readBuf, recvData);
        fftObj.fft(recvData, recvSpec);

        auto res = fftObj.findConvolutionPeak(sendSpec.frequencyDomain, recvSpec.frequencyDomain, rsltData, 40, 10, true);
        if(!res.index.isNull){
            typeof(return) dst;
            dst.snr = res.snr;
            //writefln("iteration %s: find offset: %s -> %s.", iterIdx, iterIdx*blockSize/2+res.index.get, cast(ptrdiff_t)(iterIdx*blockSize/2+res.index.get) - cast(ptrdiff_t)offsetTX);
            dst.index = Nullable!ptrdiff_t(cast(ptrdiff_t)(iterIdx*blockSize/2+res.index.get) - cast(ptrdiff_t)offsetTX);
            return dst;
        }
    }

    ConvResult!ptrdiff_t null_;
    return null_;
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
do{
    foreach(i, e; input) output[i] = complex!float(e.re, e.im);
}

struct ConvResult(Index)
{
    Nullable!Index index;
    real snr;

    real snrdB() @property { return 10*log10(snr); }
}


ConvResult!size_t findConvolutionPeak(FftObj)(
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

    //if(snrdB > 10) writeln("SNR: ", snrdB);

    if(snrdB > dBThreshold)
        return ConvResult!size_t(Nullable!size_t(maxIdx), snr);
    else
        return ConvResult!size_t(nullV, snr);
}
