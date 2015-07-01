module dffdd.filter.lms;

import std.complex;
import std.math;

final class LMSAdapter(State)
{
    import std.algorithm;

    enum size_t N = typeof(State.init.state).length;
    enum size_t P = typeof(State.init.state[0]).length;

    alias C = typeof(State.init.state[0][0]);
    alias F = typeof(State.init.power[0]);

    this(State state, real mu, size_t forgetCycle = 1024, real forgetCoeff = 0.50)
    {
        _mu = mu;
        _cnt = 0;
        _cycle = forgetCycle;
        _fcoeff = forgetCoeff;
        foreach(ref e; _subMaxPower) e = 0;
        _maxPower = state.power;
    }


    void adapt(ref State state, C error)
    {
        foreach(i; 0 .. P)
            _subMaxPower[i] = max(state.power[i], _subMaxPower[i]);

        ++_cnt;
        if(_cnt == _cycle){
            _cnt = 0;

            foreach(i; 0 .. P)
                _maxPower[i] = _maxPower[i] * _fcoeff + _subMaxPower[i] * (1 - _fcoeff);

            foreach(ref e; _subMaxPower) e = 0;
        }

        foreach(i; 0 .. N)
            foreach(p; 0 .. P)
                state.weight[i][p] += _mu * error * state.state[i][p].conj / _maxPower[p];
    }


  private:
    immutable real _mu;
    size_t _cnt;
    immutable size_t _cycle;
    immutable real _fcoeff;
    F[P] _maxPower;
    F[P] _subMaxPower;
}
