module dffdd.blockdiagram.noise;

import std.math,
       std.random;

import std.complex;

import dffdd.utils.unit;

enum real BoltzmannConst = 1.3806488e-23;

auto boxMullerNoise()
{
    return BoxMuller!Random(rndGen);
}


struct BoxMuller(RNG)
{
    this(RNG urng)
    {
        _urng = urng;
        _x = uniform01(_urng);
        _y = uniform01(_urng);
    }


    Complex!real front() const @property
    {
        return sqrt(-2 * log(_x)) * Complex!real(cos(2*PI*_y), sin(2*PI*_y));
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


    BoxMuller!RNG save() @property
    {
        typeof(return) dst;

        dst._urng = this._urng;
        dst._x = this._x;
        dst._y = this._y;
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


struct ThermalNoise
{
    this(real sampFreq, real tempK)
    {
        _rnd = boxMullerNoise();
        _gain = sqrt(noisePower(sampFreq, tempK) / 2);
        _rnd.seed(unpredictableSeed());
    }


    Complex!real front() const @property
    {
        return _rnd.front * _gain;
    }


    enum bool empty = false;


    void popFront()
    {
        _rnd.popFront();
    }


    ThermalNoise save() @property
    {
        typeof(return) dst;

        dst._rnd = this._rnd;
        dst._gain = this._gain;
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

      private:
        Rnd _rnd;
        real _front;
    }

    Uniform01Rnd res;
    res._rnd.seed(seedValue);
    res.popFront();
    return res;
}

