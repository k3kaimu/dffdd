module dffdd.filter.iir;


/**
(1 - αz^-1)^-1 というIIRフィルタを構成します
*/
struct IIRFilter1(C)
{
    this(real alpha)
    {
        _alpha = alpha;
        _last = 0;
    }


    void apply(C1 : C)(in C1[] x, C[] o)
    in{
        assert(x.length == o.length);
    }
    body{
        o[] = x[];
        foreach(ref e; o){
            auto ee = e;
            e += _last * _alpha;
            _last = ee;
        }
    }


  private:
    real _alpha;
    real _last;
}