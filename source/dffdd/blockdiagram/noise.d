module dffdd.blockdiagram.noise;

import std.math,
       std.random;

import dffdd.math;

import dffdd.utils.unit;

enum real BoltzmannConst = 1.3806488e-23;

BoxMuller!Random boxMullerNoise()
{
    Random rnd;
    rnd.seed(unpredictableSeed());

    return BoxMuller!Random(rnd);
}


struct BoxMuller(RNG, C = StdComplex!real)
{
    this(RNG urng)
    {
        _urng = urng;
        _x = uniform01(_urng);
        _y = uniform01(_urng);
    }


    C front() const @property
    {
        return cast(C)(sqrt(-2 * log(_x)) * C(cos(2*PI*_y), sin(2*PI*_y)));
    }


    enum bool empty = false;


    void popFront()
    {
        _x = uniform01(_urng);
        _y = uniform01(_urng);
    }


    void seed(uint seed)
    {
        _urng.seed(seed);
    }


    BoxMuller save() @property
    {
        typeof(return) dst = this;

        dst._urng = this._urng.save;

        return dst;
    }


  private:
    RNG _urng;
    real _x;
    real _y;
}


real noisePower(real bandWidth, real tempK)
{
    return bandWidth * tempK * BoltzmannConst;
}


struct ThermalNoise(C = StdComplex!real)
{
    this(real sampFreq, real tempK, uint s = unpredictableSeed())
    {
        _rnd = boxMullerNoise();
        _gain = sqrt(noisePower(sampFreq, tempK) / 2);
        _rnd.seed(s);
    }


    C front() const @property
    {
        return cast(C)(_rnd.front * _gain);
    }


    enum bool empty = false;


    void popFront()
    {
        _rnd.popFront();
    }


    ThermalNoise save() @property
    {
        typeof(return) dst = this;

        dst._rnd = this._rnd.save;

        return dst;
    }


  private:
    BoxMuller!Random _rnd;
    real _gain;
}


auto uniform01Range(Rnd = Random)(uint seedValue)
{
    static struct Uniform01Rnd
    {
        auto front() const @property { return _front; }
        enum bool empty = false;


        void popFront()
        {
            _front = uniform01(_rnd);
        }


      static if(isForwardRange!Rnd)
      {
        typeof(this) save() @property
        {
            typeof(return) dst = this;

            dst._rnd = this._rnd.save;

            return dst;
        }
      }


      private:
        Rnd _rnd;
        real _front;
    }

    Uniform01Rnd res;
    res._rnd.seed(seedValue);
    res.popFront();
    return res;
}

