module dffdd.blockdiagram.lna;


struct LNA(R)
{
    this(R r, Gain gain, Gain nf, SamplingFrequency sampFreq)
    {
        _r = r.connectTo!VGA(gain);
        _noise = ThermalNoise(sampFreq, 290).connectTo!VGA(sqrt(gain.gainP * (nf.gain - 1)).gain);
    }


    auto front()
    {
        return _r.front * _g1V + _noise.front;
    }


    void popFront()
    {
        _r.popFront();
        _noise.popFront();
    }


    bool empty()
    {
        return _r.empty;
    }


  static if(isForwardRange!R)
  {
    typeof(this) save() @property
    {
          typeof(return) dst = this;

          dst._r = this._r.save;
          dst._noise = this._noise.save;
          
          return dst;
    }
  }


  private:
    VGA!R _r;
    ThermalNoise _noise;
}
