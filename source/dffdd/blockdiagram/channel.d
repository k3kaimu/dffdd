module dffdd.blockdiagram.channel;

import std.complex;
import std.math;
import std.range;
import std.random;


final class ClarkesFlatFader
{
    this(size_t n, real fdts, uint rndSeed = 114514)
    {
        Random rnd;
        rnd.seed(rndSeed);
        this(n, fdts, rnd);
    }


    this(Rnd)(size_t n, real fdts, ref Rnd rnd)
    {
        immutable real th = uniform01(rnd) * 2 * PI;

        _n = n;

        foreach(i; 0 .. n)
        {
            immutable real am = uniform01(rnd) * 2 * PI;
            immutable real bm = uniform01(rnd) * 2 * PI;

            _lastICOS ~= cos(am);
            _lastISIN ~= sin(am);
            _lastQCOS ~= cos(bm);
            _lastQSIN ~= sin(bm);

            immutable real cm = cos(((2 * i + 1) * PI + th)/4/n) * fdts * 2 * PI;

            _deltaCOS ~= cos(cm);
            _deltaSIN ~= sin(cm);
        }

        this.popFront();
    }


    enum bool empty = false;


    Complex!real front() const @property
    {
        return _frnt;
    }


    void popFront()
    {
        _frnt = Complex!real(0, 0);

        foreach(i; 0 .. _n){
            immutable IC = _lastICOS[i],
                      IS = _lastISIN[i],
                      QC = _lastQCOS[i],
                      QS = _lastQSIN[i],
                      DC = _deltaCOS[i],
                      DS = _deltaSIN[i];
            
            _lastICOS[i] = IC * DC - IS * DS;
            _lastISIN[i] = IS * DC + IC * DS;

            _lastQCOS[i] = QC * DC - QS * DS;
            _lastQSIN[i] = QS * DC + QC * DS;

            _frnt += Complex!real(_lastICOS[i], _lastQSIN[i]);
        }

        _frnt /= sqrt(_n*1.0L);
    }


    ClarkesFlatFader save() const @property
    {
        ClarkesFlatFader fader = new ClarkesFlatFader(_n, 0);

        fader._lastICOS = this._lastICOS.dup;
        fader._lastISIN = this._lastISIN.dup;
        fader._lastQCOS = this._lastQCOS.dup;
        fader._lastQSIN = this._lastQSIN.dup;
        fader._deltaCOS = this._deltaCOS.dup;
        fader._deltaSIN = this._deltaSIN.dup;
        fader._frnt = this._frnt;

        return fader;
    }


  private:
    size_t _n;              // 素波数
    real[] _lastICOS,
           _lastISIN,
           _lastQCOS,
           _lastQSIN,
           _deltaCOS,
           _deltaSIN;

    Complex!real _frnt;
}


final class StopableFader(Fader)
{
    alias E = ElementType!Fader;


    this(Fader fader, const(bool)* flag)
    {
        _fader = fader;
        _flag = flag;
    }


    E front() @property
    {
        return _fader.front;
    }


    void popFront()
    {
        if(!*_flag) _fader.popFront();
    }


  static if(isInfinite!Fader)
  {
    bool empty = false;
  }
  else
  {
    bool empty() { return _fader.empty; }
  }


  static if(isForwardRange!Fader)
  {
    auto save() @property
    {
        return new StopableFader(_fader, _flag);
    }
  }

  private:
    Fader _fader;
    const(bool)* _flag;
}


final class SelectiveFadingChannel(C, Fader = InputRange!C)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(Fader[] faders)
    {
        _faders = faders;
        _coefs.length = faders.length;
        _fir = FIRFilterConverter!C(_coefs);
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        foreach(i, ref e; _coefs){
            assert(!_faders[i].empty);
            e = _faders[i].front;
            _faders[i].popFront();
        }

        _fir(input, output);
    }


    SelectiveFadingChannel dup() @property
    {
        auto dst = new SelectiveFadingChannel();
        foreach(e; this._faders)
            dst._faders ~= e.save;
        
        dst._coefs = this._coefs.dup;
        dst._fir = FIRFilterConverter!C(dst._coefs);

        return dst;
    }


  private:
    Fader[] _faders;
    C[] _coefs;
    FIRFilterConverter!C _fir;

    private this() {}
}
