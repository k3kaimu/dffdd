module dffdd.filter.rxnonlincanc;

import std.typecons;
import std.numeric : gcd;

final class TxRxSeparateCanceller(C, TxCanceller, RxCanceller)
{
    this(TxCanceller txCanceller, RxCanceller rxCanceller)
    {
        _txCanceller = txCanceller;
        _rxCanceller = rxCanceller;
    }


    size_t inputBlockLength() @property
    {
        return gcd(_txCanceller.inputBlockLength, _rxCanceller.inputBlockLength);
    }


    void apply(Flag!"learning" doLearning)(in C[] transmit, in C[] received, C[] errors)
    in{
        assert(transmit.length == received.length);
        assert(transmit.length == errors.length);
        assert(transmit.length % this.inputBlockLength == 0);
    }
    body{
        static if(doLearning) {
            _savedTrainTxSamples ~= transmit;
            _savedTrainRxSamples ~= received;

            errors[] = received[];
        } else {
            if(_savedTrainTxSamples.length != 0)
            {
                _txCanceller.apply!(Yes.learning)(transmit, received, errors);
                _txCanceller.apply!(No.learning)(transmit, received, errors);

                // _txCancellerが再現した信号を得る
                _bufferTxRegen.length = transmit.length;
                foreach(i; 0 .. transmit.length)
                    _bufferTxRegen[i] = received[i] - errors[i];
                
                _bufferRxCancelled.length = transmit.length;
                _rxCanceller.apply!(Yes.learning)(_bufferTxRegen, errors, _bufferRxCancelled);

                _savedTrainTxSamples = null;
                _savedTrainRxSamples = null;
            }

            _txCanceller.apply!(No.learning)(transmit, received, errors);
            // _txCancellerが再現した信号を得る
            _bufferTxRegen.length = transmit.length;
            foreach(i; 0 .. transmit.length)
                _bufferTxRegen[i] = received[i] - errors[i];

            _bufferRxCancelled.length = transmit.length;
            _rxCanceller.apply!(No.learning)(_bufferTxRegen, errors, _bufferRxCancelled);
            errors[] = _bufferRxCancelled[];
        }
    }


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator)
    {
        _txCanceller.preLearning(model, signalGenerator);
    }


  private:
    TxCanceller _txCanceller;
    RxCanceller _rxCanceller;
    immutable(C)[] _savedTrainTxSamples, _savedTrainRxSamples;
    C[] _bufferTxRegen;
    C[] _bufferRxCancelled;
}
