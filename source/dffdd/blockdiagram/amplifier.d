module dffdd.blockdiagram.amplifier;

import std.traits;
import std.range;
import std.math;

import dffdd.utils.unit;


struct PowerAmplifier(R)
{
    this(R r, Gain gain, Voltage iip3, Voltage iip5 = Voltage(0)/*, Voltage iip7 = Voltage(0)*/)
    {
        _r = r;

        _g1V = gain.gain;
        _g3V = iip3.V == 0 ? 0 : (gain.gain / iip3.V^^2);
        _g5V = iip5.V == 0 ? 0 : (gain.gain / iip5.V^^4);
        //_g7V = iip7.V == 0 ? 0 : (gain.gain / iip7.V^^6);
    }


    auto front()
    {
        auto x = _r.front;
        auto x1p = x.re^^2 + x.im^^2,
             x3 = x * x1p,
             x5 = x3 * x1p;
        //     x7 = x5 * x1p;

        return x * _g1V + x3 * _g3V + x5 * _g5V;// + x7 * _g7V;
    }


    void popFront()
    {
        _r.popFront();
    }


    bool empty()
    {
        return _r.empty;
    }


  static if(isForwardRange!R)
  {
    typeof(this) save() @property
    {
        typeof(this) dst;
        dst._r = this._r.save;
        dst._g1V = this._g1V;
        dst._g3V = this._g3V;
        dst._g5V = this._g5V;

        return dst;
    }
  }


  private:
    R _r;
    real _g1V;
    real _g3V;
    real _g5V;
    //real _g7V;
}


struct VGA(R)
{
    pragma(msg, R);

    this(R r, Gain gain)
    {
        _r = r;
        _gain1V = gain.gain;
    }


    auto front()
    {
        return _r.front * _gain1V;
    }


    bool empty()
    {
        return _r.empty;
    }


    void popFront()
    {
        _r.popFront();
    }


  private:
    R _r;
    real _gain1V;
}



struct PowerControlAmplifier
{
    static
    auto makeBlock(R)(R r, Voltage op, size_t avgSize = 512)
    {
        return PowerControlAmplifierImpl!R(r, op, avgSize);
    }


    static
    struct PowerControlAmplifierImpl(R)
    {
        this(R r, Voltage op, size_t avgSize)
        {
            _r = r;
            if(_r.empty){
                _empty = true;
                return;
            }

            _power = op.V^^2;
            _alpha = 1;
            _front = _r.front;
            _r.popFront();

            _empty = false;
            _cnt = 0;
            _avgSize = avgSize;
            _sumPower = 0;
            _avgCount = 0;
        }


        auto front() const @property
        {
            return _alpha * _front;
        }


        bool empty() const @property { return _empty; }


        void popFront()
        {
            if(_r.empty){
                _empty = true;
                return;
            }

            _sumPower += _front.re^^2 + _front.im^^2;
            ++_cnt;
            if(_cnt == _avgSize){
                if(_sumPower == 0)
                    _alpha = _alpha;
                else
                    _alpha = _alpha * (1 - 1.0/(1 + _avgCount))
                           + sqrt(_power / (_sumPower / _avgSize)) * (1.0/(1 + _avgCount));

                //import std.stdio;
                //writefln("alpha: %s, sumPower: %s", _alpha, _sumPower);
                _sumPower = 0;
                _cnt = 0;
                if(_avgCount < 40) ++_avgCount;
            }

            _front = _r.front;
            _r.popFront();
        }


      private:
        R _r;
        real _power;
        real _alpha;
        Unqual!(ElementType!R) _front;
        bool _empty;
        size_t _cnt;
        size_t _avgSize;
        size_t _avgCount;
        real _sumPower;
    }
}

