module dffdd.filter.ph_dcm;


import std.algorithm;
import std.complex;
import std.stdio;
import std.meta;

import dffdd.filter.state;
import dffdd.filter.ls;


auto makeParallelHammersteinWithDCM(States...)(size_t tap, size_t learningUnit, size_t maxLearningCount, States states)
{
    return new ParallelHammersteinWithDCM!States(tap, learningUnit, maxLearningCount, states);
}


final class ParallelHammersteinWithDCM(States...)
{
    enum size_t P = States.length;

    this(size_t tap, size_t learningUnit, size_t maxLearningCount, States states)
    {
        _state = states;
        _learningUnit = learningUnit;
        _maxLearningSize = maxLearningCount * learningUnit;
        _tap = tap;
        _nextLearningSize = learningUnit;
    }


    void apply(bool bLearning = true, C)(in C[] tx, in C[] rx, C[] output)
    {
      static if(bLearning)
      {
        if(!_endLearning){
            _x ~= tx;
            _y ~= rx;
            while(!_endLearning && _nextLearningSize + _tap <= _x.length)
                applyImpl();

            output[] = rx[];
            return;
        }
      }

        foreach(i; 0 .. tx.length){
            Complex!float error = rx[i];
            foreach(p, _; States){
                _state[p].update(tx[i]);
                auto ev = _state[p].error(rx[i]);
                error -= (rx[i] - ev);
            }

            output[i] = error;
        }
    }


  private:
    States _state;
    immutable size_t _learningUnit, _maxLearningSize, _tap;
    size_t _nextLearningSize;
    bool _endLearning;
    immutable(Complex!float)[] _x;
    immutable(Complex!float)[] _y;


    template AdaptorTypes(size_t n)
    {
      static if(n == States.length)
        alias AdaptorTypes = AliasSeq!();
      else
        alias AdaptorTypes = AliasSeq!(LSAdapter!(States[n]), AdaptorTypes!(n+1));
    }


    void applyImpl()
    {
        foreach(i, _; States){
            auto newstate = new FIRState!(Complex!float, true)(_tap);
            newstate.weight[] = _state[i].weight[];
            _state[i] = newstate;
        }


        auto trX = _x[_tap .. _tap + _nextLearningSize];
        auto trY = _y[_tap .. _tap + _nextLearningSize];
        AdaptorTypes!0 adaptors;
        foreach(p, _; States){
            foreach(i; 0 .. _tap)   // fill taps
                _state[p].update(_x[i]);

            adaptors[p] = lsAdapter(_state[p], trX.length);
        }

        foreach(i; 0 .. trX.length){
            Complex!float error = trY[i];
            foreach(p, _; States){
                _state[p].update(trX[i]);
                auto ev = _state[p].error(trY[i]);
                error -= (trY[i] - ev);
            }

            foreach(p, _; States)
                adaptors[p].adapt(_state[p], error);
        }

        if(_nextLearningSize == _maxLearningSize){
            _endLearning = true;
            _x = null;
            _y = null;
        }else{
            _nextLearningSize = _nextLearningSize * 2;
            if(_nextLearningSize > _maxLearningSize)
                _nextLearningSize = _maxLearningSize;
        }
    }
}
