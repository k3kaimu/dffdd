module dffdd.blockdiagram.utils;

import core.thread;

import std.range;
import std.typecons;
import std.concurrency;
import std.container;

import carbon.channel;


template connectTo(alias Block)
{
    auto connectTo(R, Params...)(R r, Params params)
    {
      static if(is(typeof((){ auto b = Block(r, params); })))                   // Block is function
        return Block(r, params);
      else static if(is(typeof((){ auto b = Block.makeBlock(r, params); })))    // Block has .makeBlock
        return Block.makeBlock(r, params);
      else static if(is(typeof((){ auto b = Block!R(r, params); })))
        return Block!R(r, params);
      else static if(is(typeof((){ auto b = new Block!R(r, params); })))
        return new Block!R(r, params);                                          // Block is class
      else
        static assert(0);
        // return Block!R(r, params);
    }
}


auto consume(R, X)(ref R range, X x)
{
    import std.range;

    return put(range, x);
}


enum bool canMakeInstrument(E, alias func) = is(typeof((){
    static struct Range
    {
        E front() { return E.init; }
        void popFront() {}
        //enum bool empty = false;
        bool empty() const { return _empty; }
        bool _empty;
    }

    Range r;
    func(r);
}));


auto makeInstrumentInNewThread(E, alias func)(Duration interval = 0.seconds)
//if(canMakeInstrument!(E, func))
{
    import core.atomic;
    import std.typecons;

    alias Ch = shared(Channel!(E));
    alias Receiver = shared(Ch.Receiver);
    alias Sender = shared(Ch.Sender);


    static final class RangeSourceFromChannel
    {
        bool empty() const pure nothrow @safe @nogc @property { return _empty; }
        E front() @property { return _front; }


        void popFront()
        {
            while(_ch.queue!E.empty){
                checkEmpty();
                if(this.empty) return;
                Thread.sleep(_interval);
            }

            _front = *_ch.pop!E;
        }


        private
        void checkEmpty()
        {
            _empty = atomicLoad(*_emptyFlag);
        }


      private:
        Receiver _ch;
        Duration _interval;
        shared(bool)* _emptyFlag;

        E _front;
        bool _empty;
    }

    static void consumer(Receiver ch, Duration interval, shared(bool)* emptyFlag)
    {
        auto r = new RangeSourceFromChannel;
        r._ch = ch;
        r._interval = interval;
        r._emptyFlag = emptyFlag;
        r.popFront();
        func(r);
    }


    static struct Result
    {
        this(Sender sender, Tid tid, shared(bool)* emptyFlag)
        {
            _ch = sender;
            _tid = tid;
            _emptyFlag = emptyFlag;
            _terminated = false;
        }


        ~this()
        {
            terminate();
        }


        void terminate()
        {
            atomicStore(*_emptyFlag, true);
            _terminated = true;
        }


        void put(E e)
        {
            if(!_terminated)
                _ch.put(e);
        }


        Tid tid() @property { return _tid; }

      private:
        Sender _ch;
        Tid _tid;
        shared(bool)* _emptyFlag;
        bool _terminated;
    }


    auto ch = channel!E;
    auto emptyFlag = new shared(bool);
    auto tid = spawn(&consumer, ch.receiver, interval, emptyFlag);

    return RefCounted!Result(ch.sender, tid, emptyFlag);
}

unittest
{
    import std.range;
    import std.algorithm;
    import core.atomic;
    import std.stdio;

    import carbon.range;

    {
        static size_t countA(int a)
        {
            static shared size_t cnt = 0;

            atomicOp!"+="(cnt, 1);
            return cnt;
        }

        auto opr = makeInstrumentInNewThread!(int, a => a.filter!"a!=0".each!countA)();
        opr.consume([0, 1, 2, 0, 0, 3]);
        Thread.sleep(1.seconds);

        assert(countA(1) == 4);
        //opr.terminate();
    }

    {
        static size_t countB(int a)
        {
            static shared size_t cnt = 0;

            atomicOp!"+="(cnt, 1);
            return cnt;
        }

        auto opr = makeInstrumentInNewThread!(int,
            a => a.filter!"a!=0".tee(makeInstrumentInNewThread!(int,
                b => b.each!countB
            )()).each!((c){})
        )();

        opr.consume([0, 1, 2, 0, 0, 3]);
        Thread.sleep(1.seconds);

        assert(countB(1) == 4);
        opr.terminate();
    }
}


static final class FiberRange(E)
{
    bool empty() const pure nothrow @safe @nogc @property { return *_emptyFlag; }
    E front() @property { return _front; }


    void popFront()
    {
        // while(_ch.empty){
            if(*_emptyFlag) return;
            Fiber.yield();
        // }

        _front = *_ch;
        // _ch.removeFront();
    }


  private:
    E* _ch;
    bool* _emptyFlag;
    E _front;
}


static struct FiberBuffer(E)
{
    this(E* ch,Fiber fiber, bool* emptyFlag)
    {
        _ch = ch;
        // _ch = new E;
        _fiber = fiber;
        _emptyFlag = emptyFlag;
    }


    ~this()
    {
        terminate();
    }


    void terminate()
    {
        *_emptyFlag = true;
        //while(_fiber.state != Fiber.State.TERM) _fiber.call();
    }


    void put(E e)
    {
        if(_fiber.state != Fiber.State.TERM)
        {
            // (*_ch).insertBack(e);
            *_ch = e;
            _fiber.call();
        }
    }

  private:
    E* _ch;
    Fiber _fiber;
    bool* _emptyFlag;
}


auto makeInstrument(E, alias func)()
{
    return .makeInstrument!E(delegate void (FiberRange!E r){ func(r); });
}



auto makeInstrument(E)(void delegate(FiberRange!E) dg)
//if(canMakeInstrument!(E, func))
{
    // auto ch = channel!E;
    // DList!E* buffer = new DList!E();
    E* buffer = new E;
    auto emptyFlag = new bool;


    void consumer()
    {
        FiberRange!E r = new FiberRange!E;
        r._ch = buffer;
        r._emptyFlag = emptyFlag;
        r._front = *buffer;
        //r.popFront();

        dg(r);
    }


    Fiber fiber = new Fiber(&consumer);

    return RefCounted!(FiberBuffer!E)(buffer, fiber, emptyFlag);
}

unittest
{
    import std.range;
    import std.algorithm;

    {
        int[] arr;
        auto opr = makeInstrument!(int, a => a.filter!"a!=0".each!((b){ arr ~= b; }));
        opr.consume([0, 1, 2, 0, 0, 3]);

        assert(arr == [1, 2, 3]);
    }

    {
        int[] arr;
        auto opr = makeInstrument!(int, a => a.filter!"a!=0".tee(makeInstrument!(int, b => b.each!((c){ arr ~= c; }))).each!((b){}));
        opr.consume([0, 1, 2, 0, 0, 3]);

        assert(arr == [1, 2, 3]);
    }
}


ForwardRange!(ElementType!R) toWrappedRange(R)(R r)
if(isForwardRange!R)
{
  static if(is(R : typeof(return)))
    return r;
  else
    return inputRangeObject(r);
}


// final class ForwardRangeObject(R) : std.range.interfaces.ForwardRange!(ElementType!R)
// {
//     this(R r)
//     {
//         _r = r;
//     }

//     ElementType!R front() @property { return _r.front; }
//     ElementType!R moveFront() { return _r.moveFront(); }
//     void popFront() { _r.popFront(); }
//     bool empty() @property { return _r.empty; }
//     ForwardRangeObject!R save() @property { return new ForwardRangeObject(_r.save); }

//   private:
//     R _r;
// }

