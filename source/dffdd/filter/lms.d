module dffdd.filter.lms;

import std.complex;
import std.math;

final class PolynomialLMS(alias polyTermABS, uint P, C = Complex!float)
if(P >= 1 && is(typeof(polyTermABS(C.init, P)) : typeof(C.init.re)))
{
    import std.algorithm;

    alias F = typeof(C.init.re);

    this(real mu, size_t forgetCycle = 1024, real iniMaxPower = 0, real forgetCoeff = 0.50)
    {
        _mu = mu;
        _cnt = 0;
        _cycle = forgetCycle;
        _fcoeff = forgetCoeff;
        _maxPower = iniMaxPower;
        foreach(p, ref ee; _maxPowerCoeff) ee = polyTermABS(_maxPower, p);
    }


    void update(C[P][] w, C[P][] x, C x00, C e)
    {
        if(x.length == 0) return;

        immutable history = x.length;
        _subMaxPower = max(x00.re^^2 + x00.im^^2, _subMaxPower);
        ++_cnt;
        if(_cnt == _cycle){
            _cnt = 0;
            _maxPower = _maxPower * _fcoeff + _subMaxPower * (1 - _fcoeff);
            _subMaxPower = 0;

            foreach(p, ref ee; _maxPowerCoeff) ee = polyTermABS(_maxPower, p);
        }

        foreach(i; 0 .. history)
            foreach(p; 0 .. P)
                w[i][p] += _mu *  e * x[i][p].conj / _maxPowerCoeff[p];
    }


  private:
    immutable real _mu;
    size_t _cnt;
    immutable size_t _cycle;
    immutable real _fcoeff;
    real _maxPower;
    real[P] _maxPowerCoeff;
    real _subMaxPower;
}
