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
import std.string;


enum size_t N_avg = 1000;

enum size_t[] delay = [585, 0, 386];

enum FREQ = 5.11E9+ 100E3;
enum LAMBDA = 2.99792458E8 / FREQ;
enum DIST = 30.0E-3;

Complex!float[] cmp0;

void main(string[] args)
{
    foreach(cnt; [100, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300]){
        Complex!float[][] arrs;
        foreach(i, fn; args[1 .. $]){
            arrs ~= fn.binaryFileReader!cfloat(1000*1000/10*cnt + delay[i]).take(N_avg*100).map!(a => complex!float(a.re, a.im)).array();
            if(arrs[$-1].length < 1000) return;
        }

        immutable N = arrs.length;

        auto mat = new Complex!float[N*N].sliced(N, N);
        foreach(i; 0 .. mat.length) foreach(j; 0 .. mat[i].length){
            mat[i, j] = arrs[i].zip(arrs[j]).map!(a => a[0]*a[1].conj).sum(complex!float(0, 0));
            mat[i, j] /= N_avg;
        }
        //writeln(mat);

        if(cmp0 is null){
            auto cmp = new Complex!float[mat.length];
            foreach(i; 0 .. mat.length){
                cmp[i] = mat[0, i].conj / mat[0, i].abs;
                //cmp[i] = 1;
            }

            cmp0 = cmp;
            continue;
        }
        //writeln(cmp);

        foreach(i; 0 .. mat.length)
            arrs[i][] /= cmp0[i] * (mat[i, i] / mat[0, 0]);

        auto mat3 = new Complex!float[N*N].sliced(N, N);
        foreach(i; 0 .. mat.length) foreach(j; 0 .. mat[i].length){
            mat3[i, j] = arrs[i].zip(arrs[j]).map!(a => a[0]*a[1].conj).sum(complex!float(0, 0));
            mat3[i, j] /= N_avg;
        }
        //writeln(mat3);
        mat = mat3;


        float[] ev = new float[N];
        heev(N, cast(cfloat*)&(mat[0, 0]), ev.ptr);
        //writefln("%s,%(%s,%),", cnt*1000*1000/10/500/1000.0, ev.map!(a => 10*log10(a)));

        File file = File("music_%s.csv".format(cnt), "w");
        foreach(th_idx; 0 .. 2001){
            auto th = (2*PI / 2000 * th_idx - PI)/2;

            Complex!float[] av = new Complex!float[N];
            foreach(i; 0 .. N)
                av[i] = std.complex.expi(DIST * sin(th) / LAMBDA * 2 * PI * i);

            writeln(av);
            writeln(ev);
            writeln(mat);

            /*
            Complex!float sum = 0;
            Complex!float dpAV = 0;
            foreach(i; 0 .. 1){
                foreach(j; 0 .. N){
                    sum += mat[j, i].conj * av[j];
                    if(i == 0)
                        dpAV += (av[j] * av[j].conj);
                }
            }

            sum *= sum.conj;

            file.writefln("%s,%s,%s,", cnt*1000*1000/10/500/1000.0, th / PI * 180, 10*log10((dpAV/sum).re));
            */
            Complex!float sum = 0;
            foreach(i; 0 .. 2){
                Complex!float subsum = 0;
                foreach(j; 0 .. N){
                    subsum += mat[j, i].conj * av[j];
                }

                sum += 1 / (subsum.re^^2 + subsum.im^^2);
            }
            file.writefln("%s,%s,%s,", cnt*1000*1000/10/500/1000.0, th / PI * 180, 10*log10(sum.re));
        }
    }
}

