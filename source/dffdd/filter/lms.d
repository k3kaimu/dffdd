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

    this(State state, real mu)
    {
        _mu = mu;
        _cnt = 0;

      static if(state.state.shape.length > 1)
      {
        _power = new F[state.power.elementsCount].sliced(state.power.shape);

        _power[] = state.power[];
      }
      else
      {
        _power = state.power;
      }
    }


    void adapt(ref State state, C error) /*pure nothrow @safe @nogc*/
    {
        static if(!is(typeof(_power) == F))
        {
          F[] ps = new F[state.power.elementsCount];
          ps[]= 0;
          foreach(j, ref p; ps)
              foreach(i; 0 .. state.numOfTaps)
                  ps[j] += state.state[i][j].sqAbs;

          immutable pmu = max(_mu, 0.01);
          ps[] += 1E-20;
          _power[] *= (1 - pmu);
          ps[] *= pmu;
          _power[] += ps[];
        }

        foreach(i; 0 .. state.numOfTaps){
          static if(is(typeof(_power) == F))
          {
            state.weight[i] += _mu * error * state.state[i].conj / state.state[i].sqAbs;
          }
          else
          {
            foreach(ref w, s, p; lockstep(state.weight[i].byElement, state.state[i].byElement, _power.byElement))
                w += _mu * error * s.conj / p;
          }
        }
    }


  private:
    immutable real _mu;
    size_t _cnt;

  static if(State.init.state.shape.length > 1)
  {
    Slice!(State.init.power.shape.length, F*) _power;
  }
  else
  {
    F _power;
  }
}


auto lmsAdapter(State)(State state, real mu)
{
    return new LMSAdapter!State(state, mu);
}
