module app;

import std.range;
import std.algorithm;
import std.stdio;
import std.experimental.ndslice;
import std.complex;
import std.numeric;

import dffdd.utils.linalg;
import dffdd.blockdiagram.inst.file;
import carbon.linear;
import dffdd.dsp.convolution;
import dffdd.utils.fft;
import std.math;


enum size_t N_avg = 1000;

enum size_t[] delay = [585, 0, 386];


void main(string[] args)
{
    /*
    immutable cnt = 1;
    Fft fft = new Fft(N_avg);
    Complex!double[][] arrs;
    foreach(i, fn; args[1 .. $]){
        arrs ~= fft.fft(fn.binaryFileReader!cfloat(1000*1000/10*cnt + delay[i]).take(N_avg).map!(a => complex!float(a.re, a.im)).array());
        if(arrs[$-1].length < 1000) return;
    }

    File("ir_01.csv", "w").writefln("%(%s,%|\n%)", fft.convolutionPower(arrs[0].frequencyDomain, (arrs[1]).frequencyDomain, null).map!"a.re^^2+a.im^^2");
    File("ir_12.csv", "w").writefln("%(%s,%|\n%)", fft.convolutionPower((arrs[1]).frequencyDomain, (arrs[2]).frequencyDomain, null).map!"a.re^^2+a.im^^2");
    File("ir_20.csv", "w").writefln("%(%s,%|\n%)", fft.convolutionPower((arrs[2]).frequencyDomain, (arrs[0]).frequencyDomain, null).map!"a.re^^2+a.im^^2");
    */

    foreach(cnt; 0 .. 400){
        Complex!float[][] arrs;
        foreach(i, fn; args[1 .. $]){
            arrs ~= fn.binaryFileReader!cfloat(1000*1000/10*cnt + delay[i]).take(N_avg).map!(a => complex!float(a.re, a.im)).array();
            if(arrs[$-1].length < 1000) return;
        }

        immutable N = arrs.length;


        auto mat = new Complex!float[N*N].sliced(N, N);
        foreach(i; 0 .. mat.length) foreach(j; 0 .. mat[i].length){
            mat[i, j] = arrs[i].zip(arrs[j]).map!(a => a[0]*a[1].conj).sum(complex!float(0, 0));
            mat[i, j] /= N_avg;
        }
        //writeln(mat);

        float[] ev = new float[N];
        heev(N, cast(cfloat*)&(mat[0, 0]), ev.ptr);
        writefln("%s,%(%s,%),", cnt*1000*1000/10/500/1000.0, ev.map!(a => 10*log10(a)));
        //writeln(mat);
    }
}
