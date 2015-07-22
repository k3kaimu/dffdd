module dffdd.filter.ls;

import std.complex;
import std.math;

import std.stdio;


final class LSAdapter(State, size_t NumOfADCBits = 12)
{
    import std.algorithm;

    enum bool usePower = false;

  private
  {
    enum size_t N = typeof(State.init.state).length;
    enum size_t M = typeof(State.init.state[0]).length;

    alias C = typeof(State.init.state[0][0]);
    alias F = typeof(State.init.state[0][0].re);

    C[] _mx, _yv;
    //immutable size_t _reComputeCycle;
    immutable size_t _L;

    size_t _fillCNT;
    //size_t _cycleCNT;
  }


    this(size_t L, /*size_t reComputeCycle*/)
    {
        _L = L;
        _mx = new C[N * M * L];
        _yv = new C[max(N*M, L)];
        _fillCNT = 0;
        //_cycleCNT = 0;

        //_reComputeCycle = reComputeCycle;
        //_cycleCNT = _reComputeCycle;
    }


    void adapt(ref State state, C error)
    {
        //if(_cycleCNT == 0){
            {
                auto p = cast(C*)state.state.ptr,
                     q = _mx.ptr;
                foreach(i; 0 .. M*N) q[_fillCNT + i*_L] = p[i];
            }
            _yv[_fillCNT] = error;
            ++_fillCNT;

            if(_fillCNT == _L){
                //writeln(_mx);
                //writeln(_mx);
                //writeln(state.weight);
                update(state);
                _fillCNT = 0;
                //_cycleCNT = _reComputeCycle;
            }
        //}else
            //--_cycleCNT;
    }


  private:
    bool update(ref State state)
    {
        if(leastSquare(cast(C*)state.weight.ptr) == 0)
            return true;
        else return false;
    }

  static
  {
    F[] //rworkSpace,
        sworkSpace;
    //C[] workSpace;
  }

    int leastSquare(C* w)
    {
        static void adjustSize(T)(ref T[] arr, size_t n) { if(arr.length < n) arr.length = n; }

        adjustSize(sworkSpace, min(_L, M*N));

        int rankN = void;

        LAPACKE_cgelss(102, cast(int)_L, cast(int)M*N, 1, _mx.ptr, cast(int)_L, _yv.ptr, cast(int)max(_L, M*N), sworkSpace.ptr, 0.001f, &rankN);

        foreach(i; 0 .. N*M) w[i] += _yv[i];
        return 0;
    }
}

alias lapack_int = int;
alias lapack_complex_float = cfloat;
extern(C) lapack_int LAPACKE_cgelss( lapack_int matrix_order, lapack_int m, lapack_int n,
                    lapack_int nrhs, lapack_complex_float* a,
                    lapack_int lda, lapack_complex_float* b,
                    lapack_int ldb, float* s, float rcond,
                    lapack_int* rank );





