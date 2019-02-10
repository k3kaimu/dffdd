module dffdd.filter.ls;

import std.complex;
import std.math;
import std.stdio;

import dffdd.utils.linalg;

import mir.ndslice;


LSAdapter!(State, NumOfADCBits) makeLSAdapter(size_t NumOfADCBits = 12, State)(State state, size_t L)
{
    return new typeof(return)(state, L);
}


final class LSAdapter(State, size_t NumOfADCBits = 12)
{
    import std.algorithm;

  private
  {
    alias C = State.StateElementType;
    alias R = typeof(C.init.re);

    Slice!(C*, 2, Contiguous) _mx;
    C[] _yv;
    immutable size_t _L;
    size_t _fillCNT;
  }


    this(State state, size_t L)
    {
        _L = L;
        _mx = new C[state.state.elementsCount * L].sliced(state.state.elementsCount, L);
        _yv = new C[max(state.state.elementsCount, L)];
        _fillCNT = 0;
    }


    void adapt()(ref State state, C error)
    {
        _mx[0 .. $, _fillCNT] = state.state.flattened;
        _yv[_fillCNT] = error;
        ++_fillCNT;

        if(_fillCNT == _L){
            update(state);
            _fillCNT = 0;
        }
    }


  private:
    bool update(ref State state)
    {
        if(auto addW = leastSquare(state.state.elementsCount)){
            state.weight[] += addW.sliced(state.weight.shape);
            return true;
        }
        else return false;
    }

  static
  {
    R[] sworkSpace;
  }

    C[] leastSquare(size_t numOfParams)
    {
        int rankN = void;
        static void adjustSize(T)(ref T[] arr, size_t n) { if(arr.length < n) arr.length = n; }

        adjustSize(sworkSpace, min(_L, numOfParams));

        static if(is(R == float))
            alias gelss = LAPACKE_cgelss;
        else
            alias gelss = LAPACKE_zgelss;

        gelss(102, cast(int)_L, cast(int)numOfParams, 1, cast(R[2]*)&(_mx[0, 0]), cast(int)_L, cast(R[2]*)_yv.ptr, cast(int)max(_L, numOfParams), sworkSpace.ptr, 0.00001f, &rankN);

        return _yv[0 .. numOfParams];
    }
}

unittest
{
    import dffdd.filter.state;

    auto fir = new MultiFIRState!(Complex!float)(4, 4);
    auto lsAdpt = makeLSAdapter(fir, 100);
}
