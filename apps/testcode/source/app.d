import std.stdio;

import dffdd.filter.orthogonalize;
import carbon.linear;
import std.numeric;

enum size_t Dim = 4;
alias C = cdouble;
alias F = double;
import std.complex;
import std.string;
import std.math;


//SMatrix!(F, Dim, Dim) diag(M)(M arr)
//{
//    typeof(return) dst;
//    foreach(i; 0 .. Dim)
//        foreach(j; 0 .. Dim){
//            dst[i, j] = i == j ? arr[i] : 0;
//        }

//    return dst;
//}
enum size_t FFTSIZE = 4*1024*1024;

void main()
{
    Fft fft = new Fft(FFTSIZE);
    auto file = File(readln().chomp(), "r");
    auto ofile = File("output.csv", "w");
    file.seek(400 * 1024 * 8);

    cfloat[] buf = new cfloat[FFTSIZE];
    buf = file.rawRead(buf);

    Complex!float[] cbuf = new Complex!float[](FFTSIZE);
    foreach(i, ref e; cbuf)
        e = Complex!float(buf[i].re, buf[i].im);

    auto res = fft.fft(cbuf).dup;
    writeln("FFT end");

    real max = 0;
    foreach(i; 0 .. FFTSIZE/2){
        res[i] = Complex!float((res[i].re^^2 + res[i].im^^2), 0);
        if(res[i].re > max) max = res[i].re;
    }

    foreach(i; 0 .. FFTSIZE/2){
        real freq = 5.0e6 / FFTSIZE * i;
        if((500e3 - 5e3) < freq && freq < (500e3 + 5e3)){
            real dB = 10*log10((res[i].re)/max);
            ofile.writefln("%7.1f,%s,%s,", freq, res[i].re, dB);
        }
    }
}
