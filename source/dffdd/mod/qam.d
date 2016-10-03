module dffdd.mod.qam;

import std.complex;
import std.mathspecial;
import std.math;
import std.numeric;


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


struct QAM
{
    alias InputElementType = ubyte;
    alias OutputElementType = Complex!float;


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


    ref Complex!float[] modulate(in ubyte[] inputs, return ref Complex!float[] outputs)
    in{
        assert(inputs.length % _k == 0);
    }
    body{
        outputs.length = inputs.length / _k;

        foreach(i; 0 .. outputs.length){
            immutable indexRe = toInt(inputs[i*_k .. i*_k + _k/2]);
            immutable indexIm = toInt(inputs[i*_k + _k/2 .. (i+1)*_k]);

            outputs[i] = Complex!float(_toSymbol[indexRe], _toSymbol[indexIm]);
        }

        return outputs;
    }


    ref ubyte[] demodulate(in Complex!float[] inputs, return ref ubyte[] outputs)
    {
        outputs.length = inputs.length * _k;

        foreach(i; 0 .. inputs.length){
            float inpRe = inputs[i].re,
                  inpIm = inputs[i].im;

            inpRe = round((inpRe * _scale - (_L - 1.0)) / -2);
            inpIm = round((inpIm * _scale - (_L - 1.0)) / -2);

            if(inpRe < 0) inpRe = 0;
            if(inpRe > _L-1) inpRe = _L-1;
            if(inpIm < 0) inpIm = 0;
            if(inpIm > _L-1) inpIm = _L-1;

            outputs[i*_k .. i*_k + _k/2] = _grayCode[cast(size_t)round(inpRe)];
            outputs[i*_k + _k/2 .. (i+1)*_k] = _grayCode[cast(size_t)round(inpIm)];
        }

        return outputs;
    }


  private:
    size_t _k, _M, _L;
    real _scale;
    ubyte[][] _grayCode;
    float[] _toSymbol;


    size_t toInt(I)(I[] bits)
    {
        size_t dst;
        foreach(e; bits)
            dst = (dst << 1) + (e ? 1 : 0);

        return dst;
    }
}

unittest
{
    import std.stdio;

    ubyte[] inps = [0, 0, 0, 0,
                    0, 0, 0, 1,
                    0, 0, 1, 0,
                    0, 0, 1, 1,
                    0, 1, 0, 0,
                    0, 1, 0, 1,
                    0, 1, 1, 0,
                    0, 1, 1, 1,
                    1, 0, 0, 0,
                    1, 0, 0, 1,
                    1, 0, 1, 0,
                    1, 0, 1, 1,
                    1, 1, 0, 0,
                    1, 1, 0, 1,
                    1, 1, 1, 0,
                    1, 1, 1, 1,];

    QAM qam = QAM(16);
    Complex!float[] dst;
    qam.modulate(inps, dst);

    ubyte[] bs;
    qam.demodulate(dst, bs);

    assert(bs == inps);
}


real ber16QAMFromSNR(real snr)
{
    return 4.0L / 4.0L * (1.0L - 1/4.0L) / 2 * erfc(sqrt(3.0L * 4 / 2 / (16 - 1) * 10.0L^^((snr-6)/10)));
}


real snrFromBer16QAM(real ber)
{
    return findRoot(delegate(real snr){ return (ber16QAMFromSNR(snr) - ber) / ber; }, 0.0L, 25.0L);
}
