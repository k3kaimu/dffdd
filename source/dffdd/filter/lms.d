module dffdd.filter.lms;


struct LMS(F = float)
{
    this(F mu, size_t lim = 1024)
    {
        _mu = mu;
        _lim = lim;
    }


    void start(R)(R w)
    {
        //if(_x2.length != w.length)
        //    _x2.length = w.length;

        //_x2[] = 0;
        _cnt = 1;
        foreach(ref e; w) e = 0;
        //_e2 = 0;
    }


    void update(R1, R2, E)(R1 w, R2 x, E e)
    {
        import std.stdio;
        //writeln(w);
        //writeln(e);

        //_e2 += e^^2/x[i];

        foreach(i; 0 .. x.length){
            //_x2[i] = (_x2[i]*_cnt + x[i]^^2)/(_cnt+1);
            w[i] += _mu * e * x[i];// / (_x2[i] + 1);
        }
        if(_cnt < _lim) ++_cnt;
    }

  private:
    //F[] _x2;
    size_t _cnt; size_t _lim;
    F _mu;
    //F _e2;
}


/*
struct RLS(F = float)
{
    this();

    void start(size_t stateSize);
    void reset();

    void update(E)(E[] w, E[] x, E[] e);
}
*/