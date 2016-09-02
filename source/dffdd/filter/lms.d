module dffdd.filter.lms;

import std.complex;
import std.math;
import std.range : lockstep;
import std.experimental.ndslice;

final class LMSAdapter(State)
{
    import std.algorithm;

    enum bool usePower = true;
    alias C = State.StateElementType;
    alias F = typeof(C.init.re);

    //enum size_t N = typeof(State.init.state).length;
    //enum size_t P = typeof(State.init.state[0]).length;

    //alias C = typeof(State.init.state[0][0]);
    //alias F = typeof(State.init.power[0]);

    this(State state, real mu, size_t forgetCycle = 1024, real forgetCoeff = 0.50)
    {
        _mu = mu;
        _cnt = 0;
        _cycle = forgetCycle;
        _fcoeff = forgetCoeff;

      static if(state.state.shape.length > 1)
      {
        _maxPower = new F[state.power.elementsCount].sliced(state.power.shape);
        _subMaxPower = new F[state.power.elementsCount].sliced(state.power.shape);

        _maxPower[] = state.power[];
        _subMaxPower[] = 0;
      }
      else
      {
        _subMaxPower = 0;
        _maxPower = state.power;
      }
    }


    void adapt(ref State state, C error) /*pure nothrow @safe @nogc*/
    {
        import std.stdio;
        //writeln(_maxPower);
        //writeln(_subMaxPower);

      static if(is(typeof(_maxPower) == F))
        _subMaxPower = max(state.power, _subMaxPower);
      else
        foreach(ref a, ref b; lockstep(state.power.byElement, _subMaxPower.byElement))
            b = max(a, b);

        ++_cnt;
        if(_cnt == _cycle){
            _cnt = 0;

          static if(is(typeof(_maxPower) == F))
          {
            _maxPower = _maxPower * _fcoeff + _subMaxPower * (1 - _fcoeff);
            _subMaxPower = 0;
          }
          else
          {
            _maxPower[] *= _fcoeff;
            _subMaxPower[] *= (1 - _fcoeff);
            _maxPower[] += _subMaxPower[];
            _subMaxPower[] = 0;
          }
        }


        foreach(i; 0 .. state.numOfTaps){
          static if(is(typeof(_maxPower) == F))
            state.weight[i] += _mu * error * state.state[i].conj / _maxPower;
          else
            foreach(ref w, s, p; lockstep(state.weight[i].byElement, state.state[i].byElement, _maxPower.byElement))
                w += _mu * error * s.conj / p;
        }
    }


  private:
    immutable real _mu;
    size_t _cnt;
    immutable size_t _cycle;
    immutable real _fcoeff;

  static if(State.init.state.shape.length > 1)
  {
    Slice!(State.init.power.shape.length, F*) _maxPower;
    Slice!(State.init.power.shape.length, F*) _subMaxPower;
  }
  else
  {
    F _maxPower;
    F _subMaxPower;
  }
}


auto lmsAdapter(State)(State state, real mu, size_t forgetCycle = 1024, real forgetCoeff = 0.50)
{
    return new LMSAdapter!State(state, mu, forgetCycle, forgetCoeff);
}
