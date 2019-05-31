module dffdd.filter.ls;

import std.complex;
import std.math;
import std.stdio;

import dffdd.utils.linalg;
import dffdd.filter.primitives : IAdapter;

import mir.ndslice;


IAdapter!(State.StateElementType) makeLSAdapter(State)(State state, size_t L)
{
    // if(state.state.elementCount * L < 1E6)
    return new LSAdapter!State(state, L);
    // else
    // return new LowMemoryLSAdapter!State(state, L);
}


final class LSAdapter(State) : IAdapter!(State.StateElementType)
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
        _mx = new C[state.state.elementCount * L].sliced(state.state.elementCount, L);
        _yv = new C[max(state.state.elementCount, L)];
        _fillCNT = 0;
    }


    void adapt(ref State state, C error)
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
        if(auto addW = leastSquare(state.state.elementCount)){
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

    auto fir = MultiFIRState!(Complex!float)(4, 4);
    auto lsAdpt = makeLSAdapter(fir, 100);
}


final class LowMemoryLSAdapter(State) : IAdapter!(State.StateElementType)
{
    private
    {
        static if(is(typeof(State.StateElementType.init.re) == float))
        {
            alias C = cfloat;
            alias R = float;
        }
        else
        {
            alias C = cdouble;
            alias R = double;
        }
    }


    this(State state, size_t L)
    {
        _Nlearn = L;
        _Nparam = state.state.elementCount;
        _fillCNT = 0;
        _partCNT = 0;

        _xx = slice!C(_Nparam, _Nparam);
        _xy = slice!C(_Nparam);
        _temp = slice!C(_Nparam);
        _partX = slice!C(1000, _Nparam);
        _partXh = slice!C(_Nparam, 1000);
        _partY = slice!C(1000);

        _partX[] = 0 + 0i;
        _xx[] = 0 + 0i;
        _xy[] = 0 + 0i;
    }


    void adapt(ref State state, State.StateElementType error)
    {
        _temp[] = state.state.flattened.map!("cast("~C.stringof~")(a.re + a.im * 1i)");
        _partX[_partCNT, 0 .. $] = _temp;
        _partY[_partCNT] = error.re + error.im * 1i;

        ++_partCNT;
        if(_partCNT == 1000) {
            flushPart();
        }

        ++_fillCNT;

        if(_fillCNT == _Nlearn){
            flushPart();
            update(state);
            _fillCNT = 0;
        }
    }


  private:
    immutable size_t _Nlearn;
    immutable size_t _Nparam;
    size_t _fillCNT;
    size_t _partCNT;
    Slice!(C*, 2, Contiguous) _xx;
    Slice!(C*, 1, Contiguous) _xy;
    Slice!(C*, 1, Contiguous) _temp, _partY;
    Slice!(C*, 2, Contiguous) _partX, _partXh;


    void update(ref State state)
    {
        import lubeck : mldivide;
        auto ans = mldivide(_xx, _xy);
        static assert(is(typeof(ans.iterator) == C*));
        assert(ans.length == _Nparam);
        auto ws = (cast(Complex!R*)ans.iterator).sliced(state.weight.shape);
        state.weight[] += ws;
    }


    void flushPart() {
        import mir.blas : gemm, gemv;
        _partXh[] = _partX.transposed;

        C one = 1.0 + 0i;
        gemm(one, _partXh, _partX, one, _xx);
        gemv(one, _partXh, _partY, one, _xy);

        _partCNT = 0;
        _partX[] = 0 + 0i;
        _partY[] = 0 + 0i;
    }
}
