module dffdd.blockdiagram.utils;

import core.thread;
import core.lifetime : forward;

import std.meta;
import std.range;
import std.typecons;
import std.concurrency;
import std.container;

import carbon.channel;


template connectTo(alias Block)
{
    auto connectTo(R, Params...)(R r, Params params)
    if(isInputRange!R)
    {
      static if(is(typeof({static assert(isConverter!Block || isOneElementConverter!Block);})))
        return RangeAdapter!(Block, R)(r, Block(params));
      else static if(is(typeof((){ auto b = Block(r, params); })))                   // Block is function
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


auto connectTo(R, Block)(R r, Block block)
{
    return RangeAdapter!(Block, R)(r, block);
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


enum isConverter(TImpl) = is(typeof((TImpl impl){
    const(TImpl.InputElementType)[] input;
    TImpl.OutputElementType[] output;
    impl(input, output);
}));


enum bool isOneElementConverter(TImpl) = is(typeof((TImpl impl){
    TImpl.InputElementType input;
    TImpl.OutputElementType output;
    impl(input, output);
}));


enum isDuplicatableConverter(Conv) 
    = (isConverter!Conv || isOneElementConverter!Conv) 
        && is(typeof((const Conv conv) {
    Conv other = conv.dup;
}));


interface IConverter(T, U, Flag!"isDuplicatable" isDup = No.isDuplicatable)
{
    alias InputElementType = T;
    alias OutputElementType = U;

    void opCall(in T[], ref U[]);

    static if(isDup)
    {
        IConverter!(T, U, isDup) dup() const;
    }
}



struct RangeAdapter(TImpl, R)
if(is(ElementType!R : const(TImpl.InputElementType)[]))
{
    this(R range, TImpl impl)
    {
        _range = range;
        _impl = impl;
    }


    this(X...)(R range, auto ref X args)
    if(X.length == 0 || is(typeof(TImpl(args))))
    {
        static if(X.length > 0) _impl = TImpl(args);
        _range = range;
    }


    // @disable
    // this(this);


    const(OutputElementType)[] front() @property
    {
        if(_frontIsComputed) return _buffer;

        auto input = _range.front;
        _impl(input, _buffer);
        _frontIsComputed = true;

        return _buffer;
    }


    void popFront()
    {
        if(!_frontIsComputed) this.front();   // 消費する
        _range.popFront();
        _frontIsComputed = false;
    }


  static if(isInfinite!R)
  {
    enum bool empty = false;
  }
  else
  {
    bool empty() @property
    {
        return _range.empty;
    }
  }


  static if(isForwardRange!R && isDuplicatableConverter!TImpl)
  {
    typeof(this) save() @property
    {
        typeof(this) dst = this;
        dst._buffer = dst._buffer.dup;
        dst._range = dst._range.save;
        dst._impl = dst._impl.dup;

        return dst;
    }
  }


  private:
    alias InputElementType = TImpl.InputElementType;
    alias OutputElementType = TImpl.OutputElementType;

    OutputElementType[] _buffer;
    bool _frontIsComputed;
    R _range;
    TImpl _impl;
}

unittest
{
    import std.typecons;
    import std.algorithm;

    static struct Impl
    {
        alias InputElementType = int;
        alias OutputElementType = int;

        this(long* ptr) { _sum = ptr; }

        void opCall(in InputElementType[] input, ref OutputElementType[] output)
        {
            import std.algorithm : sum;

            *_sum += input.sum();
            output.length = input.length;
            output[] = input[];
        }


      private:
        long* _sum;
    }


    static
    auto applyImpl(R)(R range, long* ptr)
    {
        return RangeAdapter!(Impl, R)(range, ptr);
    }


    long* ptr = new long;
    auto r = applyImpl(iota(1024).array.chunks(32), ptr).joiner;
    assert(equal(r, iota(1024)));
    assert(*ptr == 1023*512);
}


struct RangeAdapter(TImpl, R)
if(isOneElementConverter!TImpl && is(ElementType!R : TImpl.InputElementType))
{
    this(R range, TImpl impl)
    {
        _range = range;
        _impl = impl;
    }


    this(X...)(R range, auto ref X args)
    if(X.length == 0 || is(typeof(TImpl(args))))
    {
        static if(X.length > 0)
            _impl = TImpl(args);

        _range = range;
    }


    // @disable
    // this(this);


    OutputElementType front() @property
    {
        if(!_frontIsComputed){
            auto input = _range.front;
            _impl(input, _output);
            _frontIsComputed = true;
        }
        
        return _output;
    }


    void popFront()
    {
        if(!_frontIsComputed) this.front();   // 消費する
        _range.popFront();
        _frontIsComputed = false;
    }


  static if(isInfinite!R)
  {
    enum bool empty = false;
  }
  else
  {
    bool empty() @property
    {
        return _range.empty;
    }
  }


  static if(isForwardRange!R && isDuplicatableConverter!TImpl)
  {
    typeof(this) save() @property
    {
        typeof(this) dst = this;
        dst._range = dst._range.save;
        dst._impl = dst._impl.dup;

        return dst;
    }
  }


  private:
    alias InputElementType = TImpl.InputElementType;
    alias OutputElementType = TImpl.OutputElementType;

    R _range;
    TImpl _impl;
    OutputElementType _output;
    bool _frontIsComputed;
}

unittest
{
    import std.typecons;
    import std.algorithm;

    static struct Impl
    {
        alias InputElementType = int;
        alias OutputElementType = int;

        this(long* ptr) { _sum = ptr; }

        void opCall(in InputElementType input, ref OutputElementType output)
        {
            *_sum += input;
            output = input;
        }


      private:
        long* _sum;
    }

    static assert(isOneElementConverter!Impl);


    static
    auto applyImpl(R)(R range, long* ptr)
    {
        return RangeAdapter!(Impl, R)(range, ptr);
    }


    long* ptr = new long;
    auto r = applyImpl(iota(1024), ptr);
    assert(equal(r, iota(1024)));
    assert(*ptr == 1023*512);
}


struct RangeAdapter(TImpl, R)
if(!isOneElementConverter!TImpl && is(ElementType!R : TImpl.InputElementType))
{
    this(R range, TImpl impl)
    {
        _range = range;
        _impl = impl;
    }


    this(X...)(R range, auto ref X args)
    if(X.length == 0 || is(typeof(TImpl(args))))
    {
        if(X.length != 0) _impl = TImpl(args);
        _range = range;
        _buffer = new OutputElementType[1];
    }


    // @disable
    // this(this);


    const(OutputElementType)[] front() @property
    {
        if(!_frontIsComputed){
            auto input = _range.front;
            TImpl.InputElementType[1] inputBuffer;
            inputBuffer[0] = input;
            _impl(inputBuffer[0 .. 1], _buffer);
            _frontIsComputed = true;
        }
        
        return _buffer;
    }


    void popFront()
    {
        if(!_frontIsComputed) this.front();   // 消費する
        _range.popFront();
        _frontIsComputed = false;
    }


  static if(isInfinite!R)
  {
    enum bool empty = false;
  }
  else
  {
    bool empty() @property
    {
        return _range.empty;
    }
  }


  static if(isForwardRange!R && isDuplicatableConverter!TImpl)
  {
    typeof(this) save() @property
    {
        typeof(this) dst = this;
        dst._range = dst._range.save;
        dst._impl = dst._impl.dup;
        dst._buffer = dst._buffer.dup;

        return dst;
    }
  }


  private:
    alias InputElementType = TImpl.InputElementType;
    alias OutputElementType = TImpl.OutputElementType;

    R _range;
    TImpl _impl;
    OutputElementType[] _buffer;
    bool _frontIsComputed;
}

unittest
{
    import std.typecons;
    import std.algorithm;

    static struct Impl
    {
        alias InputElementType = int;
        alias OutputElementType = int;

        this(long* ptr) { _sum = ptr; }

        void opCall(in InputElementType[] input, ref OutputElementType[] output)
        {
            import std.algorithm;

            *_sum += input.sum;
            output.length = input.length;
            output[] = input[];
        }

      private:
        long* _sum;
    }


    static
    auto applyImpl(R)(R range, long* ptr)
    {
        return RangeAdapter!(Impl, R)(range, ptr);
    }


    long* ptr = new long;
    auto r = applyImpl(iota(1024), ptr).joiner;
    assert(equal(r, iota(1024)));
    import std.stdio;
    assert(*ptr == 1023*512);
}


private
bool isValidConverterConnections(Convs...)()
{
    bool ret = true;
    static foreach(i, C; Convs) {
        static if(i != 0)
            ret = ret && is(Convs[i].InputElementType == Convs[i-1].OutputElementType); 
    }

    return ret;
}

unittest
{
    static struct DummyConv(T, U) {
        alias InputElementType = T;
        alias OutputElementType = U;
    }


    static assert(isValidConverterConnections!(
        DummyConv!(int, long), DummyConv!(long, short)
    )());


    static assert(!isValidConverterConnections!(
        DummyConv!(int, short), DummyConv!(long, short)
    )());
}


private
template ConverterIOBufferTypes(Convs...)
{
    static if(Convs.length == 0)
    {
        alias ConverterIOBufferTypes = AliasSeq!();
    }
    else static if(Convs.length == 1)
    {
        alias ConverterIOBufferTypes = AliasSeq!(Convs[0].InputElementType[], Convs[0].OutputElementType[]);
    }
    else {
        alias ConverterIOBufferTypes = AliasSeq!(Convs[0].InputElementType[], ConverterIOBufferTypes!(Convs[1 .. $]));
    }
}


unittest
{
    static struct DummyConv(T, U) {
        alias InputElementType = T;
        alias OutputElementType = U;
    }


    static assert(is(ConverterIOBufferTypes!() == AliasSeq!()));
    static assert(is(ConverterIOBufferTypes!(DummyConv!(int, double)) == AliasSeq!(int[], double[])));
    static assert(is(ConverterIOBufferTypes!(DummyConv!(int, double), DummyConv!(double, short)) == AliasSeq!(int[], double[], short[])));
}


/// 複数のコンバータをつなぎ合わせます
struct ConverterGroup(Convs...)
if(isValidConverterConnections!Convs())
{
    this(Convs convs) {
        _convs = convs;
    }


    void opCall(in typeof(_buffers[0]) input, ref typeof(_buffers[$-1]) output)
    {
        if(_buffers[0].length != input.length)
            _buffers[0].length = input.length;

        _buffers[0][] = input[];

        static foreach(i, C; Convs) {
            _convs[i](_buffers[i], _buffers[i+1]);
        }

        if(output.length != _buffers[$-1].length)
            output.length = _buffers[$-1].length;
        
        output[] = _buffers[$-1][];
    }


    ref Convs[index] getConverter(size_t index)()
    {
        return _convs[index];
    }


  private:
    Convs _convs;
    ConverterIOBufferTypes!(Convs) _buffers;
}


/// ditto
ConverterGroup!Convs makeConverterGroup(Convs...)(auto ref Convs converters)
{
    return ConverterGroup!Convs(forward!converters);
}


unittest
{
    static struct DummyConv(T, U) {
        this(U g) { gain = g; }

        U gain;

        alias InputElementType = T;
        alias OutputElementType = U;

        void opCall(in T[] input, ref U[] output)
        {
            output.length = input.length;
            foreach(i, ref e; output) e = input[i];
            output[] *= gain;
        }
    }


    auto cg = makeConverterGroup(
        DummyConv!(int, int)(2),
        DummyConv!(int, float)(3),
        DummyConv!(float, double)(0.1)
    );

    int[] input = [1, 2, 3];
    double[] output;
    cg(input, output);

    import std.math;
    assert(isClose(output[0], 0.6));
    assert(isClose(output[1], 1.2));
    assert(isClose(output[2], 1.8));
}



struct ChunkedConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;

    this(size_t size)
    {
        _chunksize = size;
        _buffer.length = size;
        _buffered = 0;
    }


    void opCall(in C[] input, ref C[] output)
    {
        if(_buffered + input.length < _chunksize) {
            _buffer[_buffered .. _buffered + input.length] = input[];
            _buffered += input.length;
            output.length = 0;
            return;
        }

        immutable ntimes = (_buffered + input.length) / _chunksize;
        immutable consume = ntimes * _chunksize;
        immutable consume_input = consume - _buffered;
        output.length = ntimes * _chunksize;
        foreach(i; 0 .. _buffered)
            output[i] = _buffer[i];
        foreach(i; 0 .. consume_input)
            output[i + _buffered] = input[i];

        immutable remain = input.length - consume_input;
        _buffered = remain;
        foreach(i; 0 .. remain)
            _buffer[i] = input[consume_input + i];
    }


    ChunkedConverter!C dup() const
    {
        ChunkedConverter!C dst;

        dst._buffer = this._buffer.dup;
        dst._chunksize = this._chunksize;
        dst._buffered = this._buffered;

        return dst;
    }


  private:
    C[] _buffer;
    size_t _chunksize;
    size_t _buffered;
}

unittest
{
    auto chunker = ChunkedConverter!int(4);
    int[] input = [1, 2, 3];
    int[] output;

    chunker(input, output);
    assert(output == []);
    assert(chunker._buffered == 3);
    assert(chunker._buffer[0 .. 3] == [1, 2, 3]);
    chunker(input, output);
    assert(output == [1, 2, 3, 1]);
    assert(chunker._buffered == 2);
    assert(chunker._buffer[0 .. 2] == [2, 3]);
    chunker(input, output);
    assert(output == [2, 3, 1, 2]);
    assert(chunker._buffered == 1);
    assert(chunker._buffer[0 .. 1] == [3]);
    chunker(input, output);
    assert(output == [3, 1, 2, 3]);
    assert(chunker._buffered == 0);
    chunker(input, output);
    assert(output == []);
    assert(chunker._buffered == 3);
    assert(chunker._buffer[0 .. 3] == [1, 2, 3]);
    chunker(input ~ input, output);
    assert(output == [1, 2, 3, 1, 2, 3, 1, 2]);
    assert(chunker._buffered == 1);
    assert(chunker._buffer[0 .. 1] == [3]);
}


struct ChunkedSource(C, R)
if(is(ElementType!R : C))
{
    // alias InputElementType = C;
    alias OutputElementType = C;

    this(R range, size_t size)
    {
        _range = range;
        _chunksize = size;
    }


    void opCall(ref C[] output)
    {
        if(output.length != _chunksize)
            output.length = _chunksize;
        
        foreach(i; 0 .. _chunksize) {
            if(_range.empty) {
                output.length = 0;
                return;
            }

            output[i] = _range.front;
            _range.popFront();
        }
    }


    bool empty() { return _range.empty; }


    ChunkedSource!(C, R) dup()
    {
        ChunkedSource!(C, R) dst;

        dst._range = this._range.save;
        dst._chunksize = this._chunksize;

        return dst;
    }


  private:
    R _range;
    size_t _chunksize;
}


struct ChunkedSource(C, R)
if(is(ElementType!R : const(C)[]))
{
    // alias InputElementType = C;
    alias OutputElementType = C;


    this(R range, size_t size)
    {
        _range = range;
        _buf = ChunkedConverter!(C)(size);
    }


    void opCall(ref C[] output)
    {
        output.length = 0;
        while(output.length == 0 && !_range.empty) {
            _buf(_range.front, output);
            _range.popFront();
        }
    }


    bool empty() { return _range.empty; }


    ChunkedSource!(C, R) dup()
    {
        ChunkedSource!(C, R) dst;

        dst._range = this._range.save;
        dst._buf = this._buf.dup;

        return dst;
    }


  private:
    R _range;
    ChunkedConverter!(C) _buf;
}


auto makeChunkedSource(C, R)(R range, size_t size)
{
    return ChunkedSource!(C, R)(range, size);
}


unittest 
{
    auto r = iota(0, 10);
    auto chunked = r.makeChunkedSource!int(3);
    int[] buf;

    assert(chunked.empty == false);
    chunked(buf);
    assert(buf == [0, 1, 2]);
    assert(chunked.empty == false);
    chunked(buf);
    assert(buf == [3, 4, 5]);
    assert(chunked.empty == false);
    chunked(buf);
    assert(buf == [6, 7, 8]);
    // assert(chunked.empty == false);
    chunked(buf);
    assert(buf == []);
    assert(chunked.empty == true);
}


unittest 
{
    import std.algorithm : map;

    auto r = iota(0, 5).map!"[a*1, a*2]"();
    auto chunked = r.makeChunkedSource!int(3);
    
    int[] buf;

    assert(chunked.empty == false);
    chunked(buf);
    assert(buf == [0, 0, 1]);
    assert(chunked.empty == false);
    chunked(buf);
    assert(buf == [2, 2, 4]);
    assert(chunked.empty == false);
    chunked(buf);
    assert(buf == [3, 6, 4]);
    // assert(chunked.empty == false);
    chunked(buf);
    assert(buf == []);
    assert(chunked.empty == true);
}
