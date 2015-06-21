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

import carbon.math;
import carbon.stream;

import dffdd.dsp.convolution;


void main(string[] args)
{
    string filenameA = "received.dat", filenameB = "reference.dat", outputFilename = "conv.csv";
    real speedRatio = 1;
    size_t blockSize = 1024;

    getopt(args,
        "recData", &filenameA,
        "refData", &filenameB,
        "output", &outputFilename,
        "speedRatio", &speedRatio,
        "blockSize", &blockSize);

    Complex!double[] inputB = getRefData(filenameB, blockSize);

    outputConvolution(filenameA, inputB, blockSize, outputFilename);
}

Complex!double[] getRefData(string filename, size_t blockSize)
{
    double[2][] buf = new double[2][blockSize];
    auto file = File(filename);
    auto rs = file.rawRead(buf);
    foreach(i; 0 .. blockSize - rs.length)
        buf[rs.length+i][] = 0;
    auto ret = new typeof(return)(blockSize);
    foreach(i; 0 .. blockSize)
        ret[i] = complex!double(buf[i][0], buf[i][1]);

    return ret;
}


void outputConvolution(CPX)(string filenameA, CPX[] refData, size_t size, string ofilename)
{
    File ifile = File(filenameA);
    File file = File(ofilename, "w");
    Complex!double[] buf1 = new Complex!double[size],
                     buf2 = new Complex!double[size],
                     buf3 = new Complex!double[size];

    double[2][] rawBuf = new double[2][size];

    Fft fftObj = new Fft(size);
    fftObj.fft(refData, buf2);

    Complex!double[][] table;

    bool endSW = false;
    while(!endSW){
        auto readBuf = ifile.rawRead(rawBuf);
        immutable readLen = readBuf.length;
        foreach(i; 0 .. readLen)
            buf1[i] = complex!double(readBuf[i][0], readBuf[i][1]);

        if(readLen != size){
            endSW = true;
            foreach(ref e; buf1[readLen .. size])
                e = complex!double(0, 0);
        }

        writeln(buf1[0 .. 10]);
        writeln(buf2[0 .. 10]);
        fftObj.fft(buf1, buf3);

        fftObj.convolutionPower(buf3, buf2, buf1);
        table ~= buf1.dup;
    }

    foreach(i; 0 .. size){
        foreach(j; 0 .. table.length)
            file.writef("%s,", table[j][i].abs);

        file.writeln();
    }
}
