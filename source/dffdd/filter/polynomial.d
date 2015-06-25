module dffdd.filter.polynomial;

import std.complex;
import std.math;

import std.stdio;

import carbon.stream;

import dffdd.filter.traits;


auto polynomialFilter(alias polyTerm, uint P, C = Complex!float, U)(U updater, size_t history)
{
    return new PolynomialFilter!(polyTerm, P, C, U)(updater, history);
}


final class PolynomialFilter(alias polyTerm, uint P, C = Complex!float, U)
if(P >= 1 && is(typeof(polyTerm(C.init, P)) : C) && isUpdater!(U, P, C))
{
    this(U updater, size_t history)
    {
        _updater = updater;
        _w = new C[P][](history);
        _x = new C[P][](history);
        _size = history;

      static if(is(std.complex.Complex!F : C, F))
        C zero = C(0, 0);
      else
        C zero = 0+0i;

        foreach(i; 0 .. history) foreach(p; 0 .. P){
            _w[i][p] = zero;
            _x[i][p] = zero;
        }
    }


    void apply(in C[] tx, in C[] rx, C[] outputBuf)
    in{
        assert(tx.length == rx.length
            && tx.length == outputBuf.length);
    }
    body{
        foreach(i; 0 .. tx.length){
            // set x[0]
            foreach(p; 0 .. P)
                _x[0][p] = polyTerm(tx[i], p);

            C error = rx[i];
            foreach(j; 0 .. _size)
                foreach(p; 0 .. P)
                    error -= _w[j][p] * _x[j][p];

            _updater.update(_w, _x, tx[i], error);
            foreach_reverse(j; 1 .. _size)
                _x[j] = _x[j-1];

            outputBuf[i] = error;
        }
    }

  private:
    U _updater;
    C[P][] _w;
    C[P][] _x;
    size_t _size;
}


/+
final class PolynomialFilter(alias polyTerm, uint P, F = float, S, G, U)
if(P >= 1 && is(typeof(polyTerm(Complex!F.init, P)) : F) && isInputStream!S && isInputStream!G && isUpdater!U)
{
    import std.algorithm : min;
    import std.stdio;
    import std.complex;

    alias C = Complex!F;


    this(S tx, G rx, U updater, size_t history, size_t bufSize = 1024)
    {
        static C[] toArray(C[P][] inp) @trusted { return (cast(C*)inp.ptr)[0 .. P*inp.length]; }

        _tx = tx;
        _rx = rx;
        _updater = updater;
        _w = new F[P][](history); _wAsArray = toArray(_w);
        _x = new F[P][](history); _xAsArray = toArray(_x);
        _txBuf = new F[bufSize];
        _txBufRem = null;
        _rxBuf = new F[bufSize];
        _idx = 0;
        _size = history;

        _updater.start(_w);
        immutable zero = C(0, 0);
        foreach(ref e; _wAsArray) e = zero;
        foreach(ref e; _xAsArray) e = zero;
    }


    E[] read(E)(E[] buf)
    {
        if(_txBufRem.length != 0){
            immutable cnt = {
                size_t cnt;
                immutable ml = min(_txBufRem.length, buf.length);
                while(!_rx.empty && cnt != ml) cnt += _rx.read(_rxBuf[cnt .. ml]).length;
                return cnt;
            }();

            if(cnt == 0) return buf[0 .. 0];

            apply(_txBufRem[0 .. cnt], _rxBuf[0 .. cnt]);
            buf[0 .. cnt] = _txBufRem[0 .. cnt];

            _txBufRem = _txBufRem[(cnt == $ ? $ : cnt) .. $];

            auto next = read(buf[cnt .. $]).length;
            return buf[0 .. cnt+next];
        }else{
            auto rem = buf;
            while(rem.length){
                immutable isize = {
                    size_t cnt;
                    immutable ml = min(_txBuf.length, rem.length);
                    while(!_tx.empty && cnt != ml) cnt += _tx.read(_txBuf[cnt .. ml]).length;
                    return cnt;
                }();
                if(isize == 0) return buf[0 .. $ - rem.length];

                immutable rsize = {
                    size_t cnt;
                    while(!_tx.empty && cnt != isize) cnt += _rx.read(_rxBuf[cnt .. isize]).length;
                    return cnt;
                }();
                if(rsize == 0){
                    _txBufRem = _txBuf[0 .. isize];
                    return buf[0 .. $ - rem.length];
                }

                immutable size = min(isize, rsize);
                apply(_txBuf[0 .. size], _rxBuf[0 .. size]);
                buf[0 .. size] = _txBuf[0 .. size];
                rem = rem[size .. $];

                if(size != _txBuf.length){
                    _txBufRem = _txBuf[size .. $];
                    return buf[0 .. $ - rem.length];
                }
            }

            return buf;
        }
    }


    private
    void apply(C[] tx, C[] rx)
    in{
        assert(tx.length == rx.length);
    }
    body{
        ////foreach(i, ref e; rfr){
        ////    F[] s = _x[_idx][];
        ////    foreach(p; 0 .. P){
        ////        s[p] = e^^(p+1);
        ////    }

        ////    e = inp[i];
        ////    foreach(j; 0 .. _size){
        ////        auto jj = (_idx - j + _size) % _size;
        ////        foreach(p; 0 .. P){
        ////            e -= _w[j][p] * _x[jj][p];
        ////        }
        ////    }

        ////    // 最適化器により，パラメータを最適化する
        ////    _updater.update(Flatten!()(0, 0, 0, 0, false, _w), Flatten!()(_idx, 0, 0, 0, false, _x), e);
        ////    _idx = (_idx + 1) % _size;
        ////}
        //foreach(i, ref e; tx){
        //    F[P] s;
        //    foreach(p, ref e; 0 .. P)
        //        e = polyTerm(e, )
        //}
        foreach(i; 0 .. tx.length){
            // set x[0]
            foreach(p, ref e; 0 .. P)
                _x[0][p] = polyTerm(tx[i], p);

            C error = rx[i];
            foreach(j; 0 .. _size)
                foreach(p; 0 .. p)
                    error -= w[j][p] * x[j][p];


        }
    }


  private:
    S _tx;
    G _rx;
    U _updater;

    C[P][] _w; C[] _wAsArray;
    C[P][] _x; C[] _xAsArray;

    F[] _txBuf; F[] _txBufRem;
    F[] _rxBuf;

    size_t _idx, _size;

    static struct Flatten()
    {
        ref inout(F) front() inout { return _arr[_i][_j]; }
        ref inout(F) opIndex(size_t i) inout
        {
            size_t ii = _i + i / P;
            size_t jj = _j + i % P;
            if(jj >= P){ jj %= P; ++ii; }
            if(ii >= _arr.length) ii %= P;
            return _arr[ii][jj];
        }
        void popFront() { b = true; ++_j; if(_j == P) { _j = 0; ++_i; } if(_i == _arr.length) _i = 0; }
        bool empty() const { return b && _i == _si && _j == _sj; }
        size_t length() const { return b ? ((P * _arr.length - (_i * P + _j) + (_si * P + _sj)) % (P * _arr.length)) : (P*_arr.length); }
        alias opDollar = length;

        size_t _si, _sj, _i, _j;
        bool b;
        F[P][] _arr;
    }
}
+/


/*
auto polynomialFilter(uint P, F = float, S, G, U)(S tx, G rx, U updater, size_t history, size_t bufSize = 1024)
if(P >= 1 && isInputStream!S && isInputStream!G && isUpdater!U)
{
    static assert(isInputStream!S);
    static assert(isInputStream!G);
    static assert(isUpdater!U);
    return new PolynomialFilter!(P, F, S, G, U)(tx, rx, updater, history, bufSize);
}


final class PolynomialFilter(uint P, F = float, S, G, U)
if(P >= 1 && isInputStream!S && isInputStream!G && isUpdater!U)
{
    import std.algorithm : min;
    import std.stdio;

    this(S tx, G rx, U updater, size_t history, size_t bufSize = 1024)
    {
        _tx = tx;
        _rx = rx;
        _updater = updater;
        _w = new F[P][](history);
        _x = new F[P][](history);
        _txBuf = new F[bufSize];
        _txBufRem = null;
        _rxBuf = new F[bufSize];
        _idx = 0;
        _size = history;

        _updater.start(Flatten!()(0, 0, 0, 0, false, _w));
        foreach(i; 0 .. history)
            foreach(j; 0 .. P)
                _x[i][j] = 0;
    }


    E[] read(E)(E[] buf)
    {
        if(_txBufRem.length != 0){
            immutable cnt = {
                size_t cnt;
                immutable ml = min(_txBufRem.length, buf.length);
                while(!_rx.empty && cnt != ml) cnt += _rx.read(_rxBuf[cnt .. ml]).length;
                return cnt;
            }();

            if(cnt == 0) return buf[0 .. 0];

            apply(_txBufRem[0 .. cnt], _rxBuf[0 .. cnt]);
            buf[0 .. cnt] = _txBufRem[0 .. cnt];

            _txBufRem = _txBufRem[(cnt == $ ? $ : cnt) .. $];

            auto next = read(buf[cnt .. $]).length;
            return buf[0 .. cnt+next];
        }else{
            auto rem = buf;
            while(rem.length){
                immutable isize = {
                    size_t cnt;
                    immutable ml = min(_txBuf.length, rem.length);
                    while(!_tx.empty && cnt != ml) cnt += _tx.read(_txBuf[cnt .. ml]).length;
                    return cnt;
                }();
                if(isize == 0) return buf[0 .. $ - rem.length];

                immutable rsize = {
                    size_t cnt;
                    while(!_tx.empty && cnt != isize) cnt += _rx.read(_rxBuf[cnt .. isize]).length;
                    return cnt;
                }();
                if(rsize == 0){
                    _txBufRem = _txBuf[0 .. isize];
                    return buf[0 .. $ - rem.length];
                }

                immutable size = min(isize, rsize);
                apply(_txBuf[0 .. size], _rxBuf[0 .. size]);
                buf[0 .. size] = _txBuf[0 .. size];
                rem = rem[size .. $];

                if(size != _txBuf.length){
                    _txBufRem = _txBuf[size .. $];
                    return buf[0 .. $ - rem.length];
                }
            }

            return buf;
        }
    }


    private
    void apply(E)(E[] inp, F[] rfr)
    in{
        assert(inp.length == rfr.length);
    }
    body{
        foreach(i, ref e; rfr){
            F[] s = _x[_idx][];
            foreach(p; 0 .. P)
                s[p] = e^^(p+1);

            e = inp[i];
            foreach(j; 0 .. _size){
                auto jj = (_idx - j + _size) % _size;
                foreach(p; 0 .. P){
                    e -= _w[j][p] * _x[jj][p];
                }
            }

            // 最適化器により，パラメータを最適化する
            _updater.update(Flatten!()(0, 0, 0, 0, false, _w), Flatten!()(_idx, 0, 0, 0, false, _x), e);
            _idx = (_idx + 1) % _size;
        }
    }


  private:
    S _tx;
    G _rx;
    U _updater;

    F[P][] _w;
    F[P][] _x;

    F[] _txBuf; F[] _txBufRem;
    F[] _rxBuf;

    size_t _idx, _size;

    static struct Flatten()
    {
        ref inout(F) front() inout { return _arr[_i][_j]; }
        ref inout(F) opIndex(size_t i) inout
        {
            size_t ii = _i + i / P;
            size_t jj = _j + i % P;
            if(jj >= P){ jj %= P; ++ii; }
            if(ii >= _arr.length) ii %= P;
            return _arr[ii][jj];
        }
        void popFront() { b = true; ++_j; if(_j == P) { _j = 0; ++_i; } if(_i == _arr.length) _i = 0; }
        bool empty() const { return b && _i == _si && _j == _sj; }
        size_t length() const { return b ? ((P * _arr.length - (_i * P + _j) + (_si * P + _sj)) % (P * _arr.length)) : (P*_arr.length); }
        alias opDollar = length;

        size_t _si, _sj, _i, _j;
        bool b;
        F[P][] _arr;
    }
}

unittest
{
    import std.random;
    import std.range;
    import std.algorithm;
    import std.stdio;
    import std.math;
    import dffdd.filter.lms;

    //auto rs = repeatStream(array(map!"cast(float)"(take(rndGen, 1024))));
    real[] rnds = new real[1024];
    foreach(i, ref e; rnds) e = sin(2*PI*i/(1024*8));
    auto rs = repeatStream(rnds);
    auto amp = rs.amplifier(10);

    auto filter = polynomialFilter!(1, real)(amp, rs, LMS!(real)(0.0001), 64);

    real[] buf = new real[1024];
    foreach(i; 0 .. 1024){
        auto es = filter.read(buf);

        real sum = 0;
        foreach(e; es)
            sum += e^^2;
        writeln(sum / es.length);
    }
}
*/