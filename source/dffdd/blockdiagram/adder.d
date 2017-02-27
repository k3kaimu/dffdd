module dffdd.blockdiagram.adder;

import std.algorithm;
import std.range;


auto add(R1, R2)(R1 r1, R2 r2)
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
