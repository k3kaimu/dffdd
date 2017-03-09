module dffdd.filter.regenerator;

import std.algorithm;
import std.complex;
import std.experimental.ndslice;
import std.math;
import std.stdio;

import carbon.math : complexZero;

import dffdd.utils.fft;
import dffdd.filter.state;


final class OverlapSaveRegenerator(C)
{
    this(Slice!(2, C*) freqResponse)
    {
        _freqResponse = freqResponse;
        _fftw = makeFFTWObject!Complex(freqResponse.length!0);
        _buffer = new C[](freqResponse.length!0);
        _inputs = new C[][](freqResponse.length!1, freqResponse.length!0);
        foreach(ref es; _inputs) foreach(ref e; es) e = complexZero!C;
    }


    Slice!(2, C*) frequencyResponse() @property { return _freqResponse; }


    void apply(const(C[])[] tx, ref C[] output)
    in{
        foreach(es; tx) foreach(e; es){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    out{
        foreach(e; output){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    body{
        immutable batchLen = _freqResponse.length!0 / 2;

        output.length = tx.length;
        size_t remain = tx.length;

        auto ops = output;
        while(remain){
            auto size = min(remain, batchLen);
            applyImpl(tx[0 .. size], ops[0 .. size]);
            tx = tx[size .. $];
            ops = ops[size .. $];
            remain -= size;
        }
    }


    void applyImpl(in C[][] tx, C[] output)
    in{
        assert(tx.length <= _freqResponse.length!0);
        assert(tx.length == output.length);
        foreach(es; _inputs) foreach(e; es){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }

        foreach(es; tx) foreach(e; es){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    out{
        foreach(e; output){
            assert(!isNaN(e.re));
            assert(!isNaN(e.im));
        }
    }
    body{
        immutable size_t size = tx.length;
        output[] = complexZero!C;

        // _bufferをゼロ初期化
        foreach(ref e; _buffer) e = complexZero!C;

        foreach(p; 0 .. _freqResponse.length!1)
        {
            // _inputsをシフトして，後ろにsize分の空きを作る
            foreach(i; 0 .. _inputs[p].length - size)
                _inputs[p][i] = _inputs[p][i + size];

            // 後ろにできた空き部分にコピーする
            foreach(i, ref e; _inputs[p][$ - size .. $])
                e = tx[i][p];

            // 信号xを周波数領域に変換しXを得る
            auto ips = _fftw.inputs!float;
            auto ops = _fftw.outputs!float;
            ips[] = _inputs[p][];
            _fftw.fft!float();

            // XとHの積を計算する
            foreach(i; 0 .. _freqResponse.length)
                _buffer[i] += ops[i] * _freqResponse[i][p];
        }

        // y = IFFT{XH}の計算
        // IFFT
        _fftw.inputs!float[] = _buffer[];
        _fftw.ifft!float();
        // 結果の後ろをoutputに足す
        output[] += _fftw.outputs!float[$ - size .. $];
    }


  private:
    FFTWObject!Complex _fftw;
    Slice!(2, C*) _freqResponse;
    C[][] _inputs;
    C[] _buffer;
}


final class MultiFIRRegenerator(C)
{
    this(Slice!(2, C*) weights)
    {
        //_weight = weights;
        _state = MultiFIRState!C(weights.length!1, weights.length!0);
        _state.weights = weights;
    }


    void apply(in C[][] tx, C[] output)
    {
        foreach(i; 0 .. tx.length){
            _state.update(tx[i]);
            output[i] += _state.output;
        }
    }


  private:
    MultiFIRState!C _state;
}
