module dffdd.mod.qam;

import std.complex;
import std.mathspecial;
import std.math;
import std.numeric;

import dffdd.mod.primitives;
import dffdd.mod.bpsk;
import dffdd.mod.qpsk;
import dffdd.utils.unit;
import dffdd.math.complex : isNarrowComplex;


Bit[][] getGrayCode(size_t N) pure nothrow @safe
{
  if(N == 1)
  {
    return [[Bit(0)], [Bit(1)]];
  }
  else
  {
    Bit[][] dst = new Bit[][](2^^N, N);
    dst.length = 2^^N;
    foreach(i, ref e; dst) e[0] = i < 2^^(N-1) ? 0 : 1;

    auto m1 = getGrayCode(N-1);
    foreach(i; 1 .. N)
        foreach(j, ref e; dst) e[i] = j < 2^^(N-1) ? m1[j][i-1] : m1[$-1-(j - 2^^(N-1))][i-1];

    return dst;
  }
}


struct QAM(C)
if(isNarrowComplex!C)
{
    alias InputElementType = Bit;
    alias OutputElementType = C;
    alias R = typeof(C.init.re);


    this(uint mary)
    {
        import std.exception : assumeUnique;

        _M = mary;
        while(mary != 1){
            mary >>= 1;
            ++_k;
        }

        _L = 2^^(_k/2);

        immutable avgE = (_M - 1) * 2.0 / 3;
        _scale = sqrt(avgE);

        _grayCode = getGrayCode(_k/2).assumeUnique;

        real[] toSymbol = new real[](_grayCode.length);
        foreach(i; 0 .. toSymbol.length)
            toSymbol[toInt(_grayCode[i])] = ((_L - 1.0) - i*2) / _scale;
          
        _toSymbol = toSymbol.assumeUnique;
    }


    size_t symInputLength() const @property { return _k; }
    size_t symOutputLength() const @property { return 1; }


    ref C[] modulate(in Bit[] inputs, return ref C[] outputs) const
    in{
        assert(inputs.length % _k == 0);
    }
    do{
        outputs.length = inputs.length / _k;

        foreach(i; 0 .. outputs.length){
            immutable indexRe = toInt(inputs[i*_k .. i*_k + _k/2]);
            immutable indexIm = toInt(inputs[i*_k + _k/2 .. (i+1)*_k]);

            outputs[i] = C(_toSymbol[indexRe], _toSymbol[indexIm]);
        }

        return outputs;
    }


    ref Bit[] demodulate(in C[] inputs, return ref Bit[] outputs) const
    {
        if(outputs.length != inputs.length * _k)
            outputs.length = inputs.length * _k;

        foreach(i; 0 .. inputs.length){
            ushort sym = demodulate_symbol_impl(inputs[i]);

            foreach(j; 0 .. _k/2) {
                outputs[(i+1)*_k - 1 - j] = sym & 1;
                sym >>= 1;
            }
            foreach(j; 0 .. _k/2) {
                outputs[(i+1)*_k - _k/2 - 1 - j] = sym & 1;
                sym >>= 1;
            }
        }

        return outputs;
    }


    ref ushort[] demodulate_symbol(in C[] inputs, ref return ushort[] symbols) const
    {
        if(symbols.length != inputs.length)
            symbols.length = inputs.length;

        foreach(i; 0 .. inputs.length)
            symbols[i] = demodulate_symbol_impl(inputs[i]);

        return symbols;
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
    immutable(Bit[])[] _grayCode;
    immutable(real)[] _toSymbol;


    static size_t toInt(I)(I[] bits)
    {
        size_t dst;
        foreach(e; bits)
            dst = (dst << 1) + (e ? 1 : 0);

        return dst;
    }


    ushort demodulate_symbol_impl(C input) const
    {
        R inpRe = input.re,
          inpIm = input.im;

        inpRe = round((_L - 1) - (inpRe * _scale + (_L - 1))/2.0);
        inpIm = round((_L - 1) - (inpIm * _scale + (_L - 1))/2.0);

        if(inpRe <= 0) inpRe = 0;
        if(inpRe > _L-1) inpRe = _L-1;
        if(inpIm <= 0) inpIm = 0;
        if(inpIm > _L-1) inpIm = _L-1;

        import std.algorithm : min, max;
        immutable ptrdiff_t rIdx = min(max(0, cast(ptrdiff_t)inpRe), _L-1);
        immutable ptrdiff_t iIdx = min(max(0, cast(ptrdiff_t)inpIm), _L-1);

        auto upper = _grayCode[rIdx];
        auto lower = _grayCode[iIdx];

        ushort sym;
        foreach(i; 0 .. _k/2) {
            sym <<= 1;
            sym += upper[i].value & 1;
        }
        foreach(i; 0 .. _k/2) {
            sym <<= 1;
            sym += lower[i].value & 1;
        }

        return sym;
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
    Bit[] inbits = L1CACode(1).map!"cast(ubyte)(a > 0 ? 0 : 1)".take(1024*3).map!(a => Bit(a)).array();

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

                ushort[] refSyms;
                mod.demodulate_symbol(dst, refSyms);

                foreach(ref e; dst){
                    e += noiseGen.front * gain;
                    noiseGen.popFront();
                }

                Bit[] bs;
                ushort[] decSyms;
                mod.demodulate(dst, bs);
                mod.demodulate_symbol(dst, decSyms);

                total += inbits.length;

                size_t partcnt = 0;
                foreach(j; 0 .. inbits.length)
                    partcnt += inbits[j] == bs[j] ? 0 : 1;

                assert(partcnt == countBitError(refSyms, decSyms));

                cnt += partcnt;

                if(cnt > 10000) break;
            }

            real ber = cnt / (total * 1.0);
            assert(isClose(ber, berQAMFromSNR!scheme(dB), 0.2, 1));
        }
    }
}

unittest
{
    // constでコピーできるかのテスト
    alias C = Complex!float;
    const QAM!C qam1;
    QAM!C qam2 = qam1; // コピーできる
}
