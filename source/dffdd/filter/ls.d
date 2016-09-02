module dffdd.filter.ls;

import std.complex;
import std.experimental.ndslice;
import std.math;

import std.stdio;


LSAdapter!(State, NumOfADCBits) lsAdapter(size_t NumOfADCBits = 12, State)(State state, size_t L, size_t cnt)
{
    return new typeof(return)(state, L, cnt);
}


final class LSAdapter(State, size_t NumOfADCBits = 12)
{
    import std.algorithm;

    enum bool usePower = false;

  private
  {
    //enum size_t N = typeof(State.init.state).length;
    //enum size_t M = typeof(State.init.state[0]).length;

    //alias C = typeof(State.init.state[0][0]);
    //alias F = typeof(State.init.state[0][0].re);
    alias C = State.StateElementType;
    alias F = typeof(C.init.re);

    Slice!(2, C*) _mx;
    C[] _yv;
    //immutable size_t _reComputeCycle;
    immutable size_t _L;
    size_t _fillCNT;
    immutable size_t _updateLimit;
    size_t _remainUpdateCNT;
    //size_t _cycleCNT;
  }


    this(State state, size_t L, size_t totalUpdateCNT)
    {
        _L = L;
        _mx = new C[state.state.elementsCount * L].sliced(state.state.elementsCount, L);
        _yv = new C[max(state.state.elementsCount, L)];
        _fillCNT = 0;
        _remainUpdateCNT = totalUpdateCNT;
        _updateLimit = totalUpdateCNT;
        //_cycleCNT = 0;

        //_reComputeCycle = reComputeCycle;
        //_cycleCNT = _reComputeCycle;
    }


    void reset()
    {
        _fillCNT = 0;
        _remainUpdateCNT = _updateLimit;
    }


    void adapt(ref State state, C error)
    {
        if(_remainUpdateCNT == 0) return;

        //if(_cycleCNT == 0){
            {
                //auto q = _mx.ptr + _fillCNT;
                ////auto p = cast(C*)state.state.ptr,
                ////     q = _mx.ptr;
                ////foreach(i; 0 .. M*N) q[_fillCNT + i*_L] = p[i];
                //foreach(ref e; state.state.byElement){
                //    *q = e;
                //    q += _L;
                //}
            }
            //_mx[0 .. $, _fillCNT] = state.state[];
            _mx[0 .. $, _fillCNT] = state.state.byElement.sliced(state.state.elementsCount);
            _yv[_fillCNT] = error;
            ++_fillCNT;

            if(_fillCNT == _L){
                //writeln(_mx);
                //writeln(_mx);
                //writeln(state.weight);
                update(state);
                _fillCNT = 0;
                _remainUpdateCNT -= 1;
                //_L *= 2;
                //_cycleCNT = _reComputeCycle;
            }
        //}else
            //--_cycleCNT;
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
    F[] //rworkSpace,
        sworkSpace;
    //C[] workSpace;
  }

    C[] leastSquare(size_t numOfParams)
    {
        static void adjustSize(T)(ref T[] arr, size_t n) { if(arr.length < n) arr.length = n; }

        adjustSize(sworkSpace, min(_L, numOfParams));

        int rankN = void;

        LAPACKE_cgelss(102, cast(int)_L, cast(int)numOfParams, 1, cast(cfloat*)&(_mx[0, 0]), cast(int)_L, cast(cfloat*)_yv.ptr, cast(int)max(_L, numOfParams), sworkSpace.ptr, 0.00001f, &rankN);

        //foreach(i; 0 .. N*M) w[i] += _yv[i];
        //return 0;
        return _yv[0 .. numOfParams];
    }
}

alias lapack_int = int;
alias lapack_complex_float = cfloat;
extern(C) lapack_int LAPACKE_cgelss( lapack_int matrix_order, lapack_int m, lapack_int n,
                    lapack_int nrhs, lapack_complex_float* a,
                    lapack_int lda, lapack_complex_float* b,
                    lapack_int ldb, float* s, float rcond,
                    lapack_int* rank );





