module dffdd.blockdiagram.decimator;


import std.traits;
import std.numeric;
import std.complex;
import std.range;

import dffdd.blockdiagram.utils;
import dffdd.blockdiagram.filter;


template connectToPowerOf2Decimator(size_t N)
{
  static if(N == 1)
  {
    auto connectToPowerOf2Decimator(R)(R r, size_t len)
    {
        return r;
    }
  }
  else static if(N == 2)
  {
    auto connectToPowerOf2Decimator(R)(R r, size_t len)
    {
        return HalfDecimator!R(r, len);
    }
  }
  else
  {
    auto connectToPowerOf2Decimator(R)(R r, size_t len)
    {
        return .connectToPowerOf2Decimator!(N/2)(r.connectTo!HalfDecimator(len), len);
    }
  }
}


struct HalfDecimator(R)
{
    alias E = ElementType!R;
    alias Re = typeof(E.re);


    this(R r, size_t len)
    {
        _r = r;
        _res = new E[len/2];
        _fftRes = new Complex!Re[len];
        _buf = new Complex!Re[len];
        
        _res[] = 0+0i;
        _fftRes[] = Complex!Re(0, 0);
        _buf[] = Complex!Re(0, 0);
        
        _i = 0;
        _len = len;
        _fftObj = new Fft(len);

        update();
    }


    E front()
    {
        return _res[_i];
    }


    void popFront()
    {
        _i += 2;
        if(_i >= _len/2){
            update();
        }
    }


    bool empty() { return _bEmpty; }


    void update()
    {
        _i = 0;

        _buf[0 .. $/2] = _buf[$/2 .. $];
        foreach(i; 0 .. _len/2){
            if(_r.empty){
                _bEmpty = true;
                return;
            }

            auto e = _r.front;
            _buf[$/2 + i] = Complex!Re(e.re, e.im);
            _r.popFront();
        }

        _fftObj.fft(_buf, _fftRes);
        _fftRes[$/4 .. $-$/4] = Complex!Re(0, 0);
        _fftObj.inverseFft(_fftRes, _buf);
        foreach(i; 0 .. _len/2)
            _res[i] = _buf[i].re + _buf[i].im*1i;
    }


  private:
    R _r;
    E[] _res;
    Complex!Re[] _buf;
    Complex!Re[] _fftRes;
    Fft _fftObj;
    size_t _i;
    size_t _len;
    bool _bEmpty;
}


struct SimpleDecimator(size_t D)
{
    static
    auto makeBlock(R)(R r)
    {
        return SimpleDecimatorImpl!(R)(r);
    }


    static struct SimpleDecimatorImpl(R)
    {
        auto front() @property { return _r.front; }
        bool empty() @property { return _r.empty; }

        void popFront()
        {
            foreach(i; 0 .. D)
                if(!_r.empty) _r.popFront();
        }

      private:
        R _r;
    }
}


struct AverageDecimator(size_t D)
{
    static
    auto makeBlock(R)(R r)
    {
        return r.splitN(D).map!"a.sum()";
    }
}



struct CICFilter(size_t D, size_t N)
{
    static
    auto makeBlock(R)(R r)
    {
      static if(N == 0)
        return r.connectTo!(SimpleDecimator!D);
      else
        return r.connectTo!IIRFilter(1).connectTo!(.CICFilter!(D, N-1)).connectTo!FIRFilter(-1);
    }
}
