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

enum FREQ = 5.11E9+ 100E3;
enum LAMBDA = 2.99792458E8 / FREQ;
enum DIST = LAMBDA / 2;


void main(string[] args)
{
    foreach(cnt; 10 .. 11){
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
        //writefln("%s,%(%s,%),", cnt*1000*1000/10/500/1000.0, ev.map!(a => 10*log10(a)));

        foreach(th_idx; 0 .. 1001){
            auto th = 2*PI / 1000 * th_idx;

            Complex!float[] av = new Complex!float[N];
            foreach(i; 0 .. N)
                av[i] = std.complex.expi(DIST * sin(th) / LAMBDA * 2 * PI * i);

            writeln(av);
            writeln(ev);
            writeln(mat);

            float sum = 0;
            foreach(i; 0 .. 1){
                foreach(j; 0 .. N){
                    auto p = mat[i, j].conj * av[j];
                    sum += p.re^^2 + p.im^^2;
                }
            }

            writefln("%s,%s,%e,", cnt*1000*1000/10/500/1000.0, th, 1/sum);
        }
    }
}

