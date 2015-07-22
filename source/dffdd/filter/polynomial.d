module dffdd.filter.polynomial;

import std.complex;
import std.math;
import std.algorithm;

import std.stdio;

import carbon.stream;
import carbon.traits;

import dffdd.filter.traits;


final class PolynomialFilter(State, Adapter)
{
    alias C = typeof(State.init.state[0][0]);
    size_t cnt;

    this(State state, Adapter adapter)
    {
        _state = state;
        _adapter = adapter;
    }


    void apply(C1 : C)(in C1[] tx, in C1[] rx, C1[] outputBuf)
    in{
        assert(tx.length == rx.length
            && tx.length == outputBuf.length);
    }
    body{
        //if(cnt / (1024 * 16) % 16 != 0){
            //_state.apply(tx, rx, outputBuf);
        //}else{
            foreach(i; 0 .. tx.length){
                _state.update(tx[i]);

                C error = _state.error(rx[i]);
                outputBuf[i] = error;

                _adapter.adapt(_state, error);
            }
        //}
        cnt += tx.length;
    }


    State state() @property { return _state; }


  private:
    State _state;
    Adapter _adapter;
}
