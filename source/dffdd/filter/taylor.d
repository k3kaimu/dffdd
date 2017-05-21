/**
    paper: https://arxiv.org/pdf/1605.01345.pdf
        - title: A Linearization Technique for Self-Interference Cancellation in Full-Duplex Radios
        - authors: Arjun Nadh, Joseph Samuel, Ankit Sharma, S. Aniruddhan, Radha Krishna Ganti
*/
module dffdd.filter.taylor;


import std.typecons;

import dffdd.filter.ls;
import dffdd.filter.state;
import dffdd.filter.primitives;



final class TaylorApproximationSICanceller(C, Dist)
{
    this(Dist distorter, size_t nLearning)
    {
        _distorter = distorter;
        _differentiator = MultiFIRState!C(1 + distorter.outputDim * 3, 9);
        foreach(i; 0 .. 9){
            foreach(p; 0 .. distorter.outputDim){
                _differentiator.weight[i, p*3] = [0,0,0,0,1,0,0,0,0][i];                  /* x */
                _differentiator.weight[i, p*3+1] = [3,-32,168,-672,0,672,-168,32,-3][i];    /* x' */
                _differentiator.weight[i, p*3+2] = [1,4,4,-4,10,-4,4,4,1][i];               /* x'' */
                _differentiator.weight[i, p*3+1] /= 840;
                _differentiator.weight[i, p*3+2] /= 64;
            }

            _differentiator.weight[i, $-1] = [0,0,0,0,1,0,0,0,0][i];                  /* y */
        }

        _cancellerState = MultiFIRState!C(distorter.outputDim * 3, 1);    /* P * count(x, x', x'') */
        _adapter = makeLSAdapter(_cancellerState, nLearning);
        _buffer = new C[][](this.inputBlockLength, distorter.outputDim);
        _ncnt = 0;
    }


    size_t inputBlockLength() @property
    {
        return _distorter.inputBlockLength;
    }


    Dist distorter() @property { return _distorter; }


    size_t delay() const @property { return 4; }


    void apply(Flag!"learning" doLearning)(in C[] input, in C[] desired, C[] errors)
    {
        immutable size_t blk = this.inputBlockLength;

        C[] diffs = new C[3 * _distorter.outputDim + 1];
        foreach(i; 0 .. input.length / blk){
            auto ips = input[i*blk .. (i+1) * blk];
            auto dss = desired[i*blk .. (i+1) * blk];
            auto ers = errors[i*blk .. (i+1) * blk];

            _distorter(ips, _buffer);
            foreach(j, e; _buffer){
                foreach(k, ee; e) foreach(n; 0 .. 3)
                    diffs[k*3 + n] = ee;
                diffs[$-1] = dss[j];
                _differentiator.update(diffs);
                foreach(k; 0 .. diffs.length)
                    diffs[k] = _differentiator.subFIRState(k).output;

                _cancellerState.update(diffs[0 .. $-1]);
                ers[j] = diffs[$-1] - _cancellerState.output;

                static if(doLearning)
                {
                    if(_ncnt >= 4)
                        _adapter.adapt(_cancellerState, ers[j]);
                    else
                        ++_ncnt;
                }
            }
        }
    }


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator)
    // if(isModelParameterSet!M)
    {
        static if(is(typeof((Dist dist, C[] txs){ dist.learn(txs); })))
        {
            auto sig = signalGenerator(model);
            auto buf = new C[](model.orthogonalizer.numOfTrainingSymbols);
            sig.fillBuffer!(["txBaseband"])(buf);
            _distorter.learn(buf);
        }
    }


  private:
    Dist _distorter;
    MultiFIRState!C _differentiator;
    MultiFIRState!C _cancellerState;
    LSAdapter!(MultiFIRState!C) _adapter;
    C[][] _buffer;
    size_t _ncnt;
}
