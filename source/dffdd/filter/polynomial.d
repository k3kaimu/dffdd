module dffdd.filter.polynomial;

import carbon.stream;

import dffdd.filter.traits;

auto polynomialFilter(uint P, F = float, S, G, U)(S stream, G reference, U updater, size_t history, size_t bufSize = 1024)
if(P >= 1 && isInputStream!S && isInputStream!G && isUpdater!U)
{
    static assert(isInputStream!S);
    static assert(isInputStream!G);
    static assert(isUpdater!U);
    return new PolynomialFilter!(P, F, S, G, U)(stream, reference, updater, history, bufSize);
}


final class PolynomialFilter(uint P, F = float, S, G, U)
if(P >= 1 && isInputStream!S && isInputStream!G && isUpdater!U)
{
    import std.algorithm : min;
    import std.stdio;

    this(S stream, G reference, U updater, size_t history, size_t bufSize = 1024)
    {
        _inp = stream;
        _rfr = reference;
        _updater = updater;
        _w = new F[P][](history);
        _x = new F[P][](history);
        _inpBuf = new F[bufSize];
        _inpBufRem = null;
        _rfrBuf = new F[bufSize];
        _idx = 0;
        _size = history;

        _updater.start(Flatten!()(0, 0, 0, 0, false, _w));
        foreach(i; 0 .. history)
            foreach(j; 0 .. P)
                _x[i][j] = 0;
    }


    E[] read(E)(E[] buf)
    {
        if(_inpBufRem.length != 0){
            immutable cnt = {
                size_t cnt;
                immutable ml = min(_inpBufRem.length, buf.length);
                while(!_rfr.empty && cnt != ml) cnt += _rfr.read(_rfrBuf[cnt .. ml]).length;
                return cnt;
            }();

            if(cnt == 0) return buf[0 .. 0];

            apply(_inpBufRem[0 .. cnt], _rfrBuf[0 .. cnt]);
            buf[0 .. cnt] = _inpBufRem[0 .. cnt];

            _inpBufRem = _inpBufRem[(cnt == $ ? $ : cnt) .. $];

            auto next = read(buf[cnt .. $]).length;
            return buf[0 .. cnt+next];
        }else{
            auto rem = buf;
            while(rem.length){
                immutable isize = {
                    size_t cnt;
                    immutable ml = min(_inpBuf.length, rem.length);
                    while(!_inp.empty && cnt != ml) cnt += _inp.read(_inpBuf[cnt .. ml]).length;
                    return cnt;
                }();
                if(isize == 0) return buf[0 .. $ - rem.length];

                immutable rsize = {
                    size_t cnt;
                    while(!_inp.empty && cnt != isize) cnt += _rfr.read(_rfrBuf[cnt .. isize]).length;
                    return cnt;
                }();
                if(rsize == 0){
                    _inpBufRem = _inpBuf[0 .. isize];
                    return buf[0 .. $ - rem.length];
                }

                immutable size = min(isize, rsize);
                apply(_inpBuf[0 .. size], _rfrBuf[0 .. size]);
                buf[0 .. size] = _inpBuf[0 .. size];
                rem = rem[size .. $];

                if(size != _inpBuf.length){
                    _inpBufRem = _inpBuf[size .. $];
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
    S _inp;
    G _rfr;
    U _updater;

    F[P][] _w;
    F[P][] _x;

    F[] _inpBuf; F[] _inpBufRem;
    F[] _rfrBuf;

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
