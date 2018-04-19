module dffdd.mod.qam;

import std.complex;
import std.mathspecial;
import std.math;
import std.numeric;

import dffdd.mod.bpsk;
import dffdd.mod.qpsk;
import dffdd.utils.unit;


ubyte[][] getGrayCode(size_t N)
{
  if(N == 1)
  {
    return cast(ubyte[][])[[false], [true]];
  }
  else
  {
    ubyte[][] dst = new ubyte[][](2^^N, N);
    dst.length = 2^^N;
    foreach(i, ref e; dst) e[0] = i < 2^^(N-1) ? 0 : 1;

    auto m1 = getGrayCode(N-1);
    foreach(i; 1 .. N)
        foreach(j, ref e; dst) e[i] = j < 2^^(N-1) ? m1[j][i-1] : m1[$-1-(j - 2^^(N-1))][i-1];

    return dst;
  }
}


struct QAM(C)
{
    alias InputElementType = ubyte;
    alias OutputElementType = C;


    this(uint mary)
    {
        _M = mary;
        while(mary != 1){
            mary >>= 1;
            ++_k;
        }

        _L = 2^^(_k/2);

        immutable avgE = (_M - 1) * 2.0 / 3;
        _scale = sqrt(avgE);

        _grayCode = getGrayCode(_k/2);

        _toSymbol.length = _grayCode.length;
        foreach(i; 0 .. _toSymbol.length)
            _toSymbol[toInt(_grayCode[i])] = ((_L - 1.0) - i*2) / _scale;
    }


    size_t symInputLength() const @property { return _k; }
    size_t symOutputLength() const @property { return 1; }


    ref C[] modulate(in ubyte[] inputs, return ref C[] outputs)
    in{
        assert(inputs.length % _k == 0);
    }
    body{
        outputs.length = inputs.length / _k;

        foreach(i; 0 .. outputs.length){
            immutable indexRe = toInt(inputs[i*_k .. i*_k + _k/2]);
            immutable indexIm = toInt(inputs[i*_k + _k/2 .. (i+1)*_k]);

            outputs[i] = C(_toSymbol[indexRe], _toSymbol[indexIm]);
        }

        return outputs;
    }


    ref ubyte[] demodulate(in C[] inputs, return ref ubyte[] outputs)
    {
        outputs.length = inputs.length * _k;

        foreach(i; 0 .. inputs.length){
            R inpRe = inputs[i].re,
              inpIm = inputs[i].im;

            inpRe = round((_L - 1) - (inpRe * _scale + (_L - 1))/2.0);
            inpIm = round((_L - 1) - (inpIm * _scale + (_L - 1))/2.0);

            if(inpRe <= 0) inpRe = 0;
            if(inpRe > _L-1) inpRe = _L-1;
            if(inpIm <= 0) inpIm = 0;
            if(inpIm > _L-1) inpIm = _L-1;

            //import std.stdio;
            //writeln(inpRe, ", ", inpIm);

            outputs[i*_k .. i*_k + _k/2] = _grayCode[cast(size_t)inpRe];
            outputs[i*_k + _k/2 .. (i+1)*_k] = _grayCode[cast(size_t)inpIm];
        }

        return outputs;
    }


    Voltage outputVoltage() const @property
    {
        real p = 0;
        foreach(e; _toSymbol)
            p += e.sqAbs();

        return Voltage(sqrt(p * 2));
    }


  private:
    size_t _k, _M, _L;
    real _scale;
    ubyte[][] _grayCode;
    real[] _toSymbol;


    size_t toInt(I)(I[] bits)
    {
        size_t dst;
        foreach(e; bits)
            dst = (dst << 1) + (e ? 1 : 0);

        return dst;
    }
}


real berQAMFromSNR(real snr, size_t M_ = 16)
{
    size_t k;
    while(M_ != 1){
        M_ >>= 1;
        ++k;
    }

    //return 4.0L / 4.0L * (1.0L - 1/4.0L) / 2 * erfc(sqrt(3.0L * 4 / 2 / (16 - 1) * 10.0L^^((snr)/10) / 4));
    immutable ebno = 10.0^^(snr/10) / k,
              M = 2^^k,
              x = sqrt(3*k*ebno/(M - 1)),
              pb = (4.0/k) * (1.0 - 1.0/sqrt(M*1.0)) * (1.0/2) * erfc(x / sqrt(2.0));

    return pb;
}


real snrFromBerQAM(real ber, size_t M = 16)
{
    return findRoot(delegate(real snr){ return (berQAMFromSNR(snr, M) - ber) / ber; }, 0.0L, 25.0L);
}


auto makeQAM(string scheme, C = Complex!float, T...)(T M)
{
  static if(scheme == "BPSK")
    return BPSK!C();
  else static if(scheme == "QPSK")
    return QPSK!C();
  else static if(scheme == "QAM")
    return QAM!C(M);
  else static if(scheme == "16QAM")
    return QAM!C(16);
  else static if(scheme == "64QAM")
    return QAM!C(64);
  else static if(scheme == "256QAM")
    return QAM!C(256);
  else
    static assert("Unsupported modulation: " ~ scheme);
}


real berQAMFromSNR(string scheme, T...)(real snr, T M)
{
  static if(scheme == "BPSK")
    return berBPSKFromSNR(snr);
  else static if(scheme == "QPSK")
    return berQPSKFromSNR(snr);
  else static if(scheme == "QAM")
    return berQAMFromSNR(snr, M);
  else static if(scheme == "16QAM")
    return berQAMFromSNR(snr, 16);
  else static if(scheme == "64QAM")
    return berQAMFromSNR(snr, 64);
  else static if(scheme == "256QAM")
    return berQAMFromSNR(snr, 256);
  else
    static assert("Unsupported modulation: " ~ scheme);
}


// BERのテスト
unittest
{
    import std.stdio;
    import std.random;
    import std.algorithm;
    import dffdd.gps.code;
    import dffdd.blockdiagram.noise;
    import dffdd.mod.bpsk;
    import std.meta;
    import std.range, std.array;

    Complex!float[] dst;

    auto noiseGen = boxMullerNoise();
    ubyte[] inbits = L1CACode(1).map!"cast(ubyte)(a > 0 ? 0 : 1)".take(1024*3).array();

    foreach(scheme; AliasSeq!("BPSK", "QPSK", "16QAM", "64QAM", "256QAM"))
    {
        int dBStart = -2,
             dBEnd = 8;

        if(scheme == "BPSK"){
            dBStart = -5;
            dBEnd = 5;
        }
        else if(scheme == "16QAM"){
            dBStart = 5;
            dBEnd = 15;
        }
        else if(scheme == "64QAM"){
            dBStart = 10;
            dBEnd = 20;
        }
        else if(scheme == "256QAM"){
            dBStart = 20;
            dBEnd = 30;
        }

        auto mod = makeQAM!scheme();

        foreach(dB; dBStart .. dBEnd+1)
        {
            immutable gain = sqrt(10.0^^(-dB/10.0) / 2);
            ulong cnt;
            ulong total;
            foreach(i; 0 .. 5_000)
            {
                mod.modulate(inbits, dst);

                foreach(ref e; dst){
                    e += noiseGen.front * gain;
                    noiseGen.popFront();
                }

                ubyte[] bs;
                mod.demodulate(dst, bs);

                total += inbits.length;
                foreach(j; 0 .. inbits.length)
                    cnt += inbits[j] == bs[j] ? 0 : 1;

                if(cnt > 10000) break;
            }

            real ber = cnt / (total * 1.0);
            assert(approxEqual(ber, berQAMFromSNR!scheme(dB), 0.2, 1));
        }
    }
}
