module dffdd.blockdiagram.adder;

import std.algorithm;
import std.range;

import rx;
import dffdd.blockdiagram.utils;

auto add(R1, R2)(R1 r1, R2 r2)
if(isInputRange!R1 && isInputRange!R2)
{
    static struct AdderResult
    {
        auto front() { return _r1.front + _r2.front; }
        void popFront()
        {
            _r1.popFront();
            _r2.popFront();
        }


        bool empty()
        {
            import std.stdio;
            return _r1.empty || _r2.empty;
        }

      static if(isForwardRange!R1 && isForwardRange!R2)
      {
        AdderResult save() @property 
        {
            AdderResult dst = this;

            dst._r1 = this._r1.save;
            dst._r2 = this._r2.save;

            return dst;
        }
      }

      private:
        R1 _r1;
        R2 _r2;
    }

    static if(is(typeof(r1 is null)))
        assert(r1 !is null);

    static if(is(typeof(r2 is null)))
        assert(r2 !is null);

    return AdderResult(r1, r2);
}


struct AdderConverterImpl(C, R)
if(isInfinite!R)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(R range)
    {
        _range = range;
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        output = input + _range.front;
        _range.popFront();
    }


  static if(isForwardRange!R)
  {
    auto dup()
    {
        return typeof(this)(_range.save);
    }
  }


  private:
    R _range;
}


alias AdderConverter(C, R) = ArrayConverterOf!(AdderConverterImpl!(C, R));


AdderConverter!(C, R) makeAdder(C, R)(R r)
if(isInfinite!R)
{
    return AdderConverter!(C, R)(r);
}


AdderConverter!(ElementType!R, R) makeAdder(R)(R r)
{
    return makeAdder!(ElementType!R, R)(r);
}


unittest
{
    import dffdd.blockdiagram.utils;

    auto r1 = iota(10);
    auto r2 = repeat(5);

    auto r3 = r1.connectTo(makeAdder(r2));
    assert(equal(iota(5, 15), r3));
}
