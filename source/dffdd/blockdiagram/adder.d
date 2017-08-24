module dffdd.blockdiagram.adder;

import std.algorithm;
import std.range;
import std.meta;
import std.variant;
import carbon.templates;

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


final class ZipAddBlock(C, size_t N) : IMNSPBlock!(TGroup!(Repeat!(N, C)), TGroup!C, true)
{
    alias ElementType = C;


    this() {
        _subject = new SubjectObject!(const(C)[])();

        foreach(i; ToTuple!(TRIota!(0, N)))
            _sinks[i] = observerObject!(const(C)[])(ObserverResult!i(this));
    }


    Variant opIndex(string str)
    {
        assert(0);
    }


    void opIndexAssignImpl(Variant v, string str)
    {
        assert(0);
    }


    bool hasSetter(string) { return false; }
    bool hasGetter(string) { return false; }


    auto subscribe(TObserver)(TObserver observer)
    {
        return doSubscribe(_subject, observer);
    }


    Sinks!(Repeat!(N, C)) sinks() @property
    {
        return _sinks;
    }


    Sources!(C) sources() @property
    {
        typeof(return) dst;
        dst[0] = this._subject;
        return dst;
    }


    ZipAddBlock dup() @property
    {
        ZipAddBlock dst = new ZipAddBlock();
        foreach(i; 0 .. N) dst._queues[i] = this._queues[i].dup;

        return dst;
    }


  private:
    C[][N] _queues;
    SubjectObject!(const(C)[]) _subject;
    Sinks!(Repeat!(N, C)) _sinks;
    C[] _buffer;


    static struct ObserverResult(size_t i)
    {
        ZipAddBlock _this;

        void put(in C[] input)
        {
            _this._queues[i] ~= input;

            size_t s = size_t.max;
            foreach(j; 0 .. N)
                s = min(s, _this._queues[j].length);

            _this._buffer.length = s;
            _this._buffer[] = C(0);

            foreach(j; 0 .. N){
                foreach(k; 0 .. s)
                    _this._buffer[k] += _this._queues[j][k];

                _this._queues[j].popFrontN(s);
            }

            _this._subject.put(_this._buffer);
        }
    }
}


auto makeZipAddBlock(C, size_t N)()
{
    return new ZipAddBlock!(C, N)();
}


unittest
{
    auto block = makeZipAddBlock!(int, 3);
    
    auto s0 = new SubjectObject!(const(int)[]);
    auto s1 = new SubjectObject!(const(int)[]);
    auto s2 = new SubjectObject!(const(int)[]);

    s0.doSubscribe(block.sinks[0]);
    s1.doSubscribe(block.sinks[1]);
    s2.doSubscribe(block.sinks[2]);

    int[] arr;
    block.doSubscribe((const(int)[] input){ arr ~= input; });

    assert(arr.empty);
    .put(s0, 1);
    assert(arr.empty);
    assert(block._queues[0] == [1]);
    assert(block._queues[1] == []);
    assert(block._queues[2] == []);
    .put(s1, 1);
    assert(arr.empty);
    assert(block._queues[0] == [1]);
    assert(block._queues[1] == [1]);
    assert(block._queues[2] == []);
    .put(s2, 1);
    assert(arr == [3]);
    assert(block._queues[0] == []);
    assert(block._queues[1] == []);
    assert(block._queues[2] == []);

    arr = arr[1 .. $];

    .put(s0, [1, 2, 3]);
    .put(s1, [4, 5]);
    .put(s2, [6, 7, 8, 9]);

    assert(arr == [11, 14]);
    assert(block._queues[0] == [3]);
    assert(block._queues[1] == []);
    assert(block._queues[2] == [8, 9]);

    arr = arr[2 .. $];

    .put(s1, [1, 2]);
    .put(s0, [1]);

    assert(arr == [12, 12]);

    foreach(i; 0 .. 3)
        assert(block._queues[i].length == 0);
}
