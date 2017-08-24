module dffdd.blockdiagram.utils;

import core.thread;

import std.range;
import std.typecons;
import std.concurrency;
import std.container;
import std.meta;

import carbon.channel;

import rx;


template connectTo_old(alias Block)
{
    auto connectTo_old(R, Params...)(R r, Params params)
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

template connectTo(alias Block)
{
    auto connectTo(R, Params...)(R r, Params params)
    if(is(typeof((){ connectTo_old!Block(r, params); })))
    {
        return connectTo_old!Block(r, params);
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


    ~this() @nogc nothrow
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


enum bool isConverter(TImpl) = is(typeof((TImpl impl){
    const(TImpl.InputElementType)[] input;
    TImpl.OutputElementType[] output;
    impl(input, output);
}));


enum bool isOneElementConverter(TImpl) = is(typeof((TImpl impl){
    TImpl.InputElementType input;
    TImpl.OutputElementType output;
    impl(input, output);
}));


enum bool isHighOrderConverter(TImpl) = is(typeof((TImpl impl){
    const(TImpl.InputElementType)[] input;
    TImpl.OutputElementType[] output;
    void delegate() flusher;
    impl(input, output, flusher);
}));


enum bool isHighOrderOneElementConverter(TImpl) = is(typeof((TImpl impl){
    const(TImpl.InputElementType) input;
    TImpl.OutputElementType output;
    void delegate() flusher;
    impl(input, output, flusher);
}));


interface IBlock(bool isDuplicable = false)
{
  static if(isDuplicable)
  {
    IBlock dup() @property;
  }

    void opIndexAssignImpl(Variant value, string key);
    Variant opIndex(string key);
    bool hasGetter(string key) @property;
    bool hasSetter(string key) @property;

    void opIndexAssign(T)(T value, string key)
    {
        static if(is(T : Variant))
            this.opIndexAssignImpl(value, key);
        else
            this.opIndexAssignImpl(Variant(value), key);
    }
}


interface ISPBlock(I, O, bool isDuplicable = false) : IBlock!(isDuplicable), IMNSPBlock!(TGroup!(I), TGroup!(O), isDuplicable)
{
    alias ElementType = const(O)[];
    void put(const(I)[] input);
    Observable!(const(O)[]) observableObject() @property;
}


struct TGroup(T...)
{
    alias Types = T;
}


alias Sinks(T...) = Tuple!(_ApplyConstArrayOf!(Observer, T));
alias Sources(T...) = Tuple!(_ApplyConstArrayOf!(Observable, T));

template _ApplyConstArrayOf(alias f, T...)
{
    alias _ApplyConstArrayOf = _ApplyConstArrayOfImpl!0;

    template _ApplyConstArrayOfImpl(size_t i)
    {
        static if(i == T.length)
            alias _ApplyConstArrayOfImpl = AliasSeq!();
        else{
            alias U = T[i];
            alias _ApplyConstArrayOfImpl = AliasSeq!(f!(const(U)[]), _ApplyConstArrayOfImpl!(i+1));
        }
    }
}


interface IMNSPBlock(Is, Os, bool isDuplicable = false) : IBlock!(isDuplicable)
if(is(Is : TGroup!X, X...) && is(Os : TGroup!Y, Y...))
{
    Sinks!(Is.Types) sinks() @property;
    Sources!(Os.Types) sources() @property;
}


class RxBlock(TImpl) :
    Observer!(const(TImpl.OutputElementType)[]),
    ISPBlock!(TImpl.InputElementType, TImpl.OutputElementType, isDuplicable!(TImpl))
if(isConverter!TImpl || isHighOrderConverter!TImpl)
{
    import std.algorithm : min;
    import std.functional : forward;
    import rx;
    import std.conv : to;
    import std.variant : Variant;

    alias InputElementType = TImpl.InputElementType;
    alias OutputElementType = TImpl.OutputElementType;

    alias ElementType =const(OutputElementType)[];


    this(TImpl impl)
    {
        _impl = impl;
        this();
    }


    this(X...)(auto ref X args)
    if(is(typeof(TImpl(args))))
    {
        this(TImpl(args));
    }


    // @disable
    // this(this);


    bool hasGetter(string key)
    {
        foreach(mname; __traits(allMembers, TImpl))
        {
            static if(
                is(typeof((){ typeof(mixin("_impl." ~ mname)) x; }))
            )
            {
                if(mname == key)
                    return true;
            }
        }

        return false;
    }


    Variant opIndex(string key)
    {
        foreach(mname; __traits(allMembers, TImpl))
        {
            static if(
                is(typeof((){ typeof(mixin("_impl." ~ mname)) x; }))
            )
            {
                if(mname == key)
                    return Variant(mixin("_impl." ~ mname));
            }
        }

        {
            import core.exception;
            onRangeError();

            Variant null_;
            return null_;
        }
    }


    bool hasSetter(string key)
    {
        foreach(mname; __traits(allMembers, TImpl))
        {
            static if(
                is(typeof((){ mixin("_impl." ~ mname); }))
              &&
                is(typeof((typeof(mixin("_impl." ~ mname)) v, TImpl impl){ mixin("impl." ~ mname ~ "=v;"); }))
              &&
                is(typeof((typeof(mixin("_impl." ~ mname)) v1, T v2){ auto x = to!(typeof(v1))(v2); }))
            )
            {
                if(mname == key)
                    return true;
            }
        }

        return false;
    }


    void opIndexAssignImpl(Variant value, string key)
    {
        foreach(mname; __traits(allMembers, TImpl))
        {
            static if(
                is(typeof((){ mixin("_impl." ~ mname); }))
              &&
                is(typeof((typeof(mixin("_impl." ~ mname)) v, TImpl impl){ mixin("impl." ~ mname ~ "=v;"); }))
              &&
                is(typeof((typeof(mixin("_impl." ~ mname)) v1, T v2){ auto x = to!(typeof(v1))(v2); }))
            )
            {
                if(mname == key){
                    alias T = typeof(mixin("_impl." ~ mname));
                    mixin("_impl" ~ mname ~ " = to!T(value)");
                    return;
                }
            }
        }

        {
            import core.exception;
            onRangeError();
        }
    }


    void opIndexAssign(T)(T value, string key)
    {
        static if(is(T : Variant))
            this.opIndexAssignImpl(value, key);
        else
            this.opIndexAssignImpl(Variant(value), key);
    }


    void put(const(InputElementType)[] input)
    {
        static if(isHighOrderConverter!TImpl)
        {
            _impl(input, _buffer, (){ .put(_subject, _buffer); });
        }
        else
        {
            _impl(input, _buffer);
            .put(_subject, _buffer);
        }
    }


    void completed() {}
    void failure(Exception ex) {}


    Observable!(const(OutputElementType)[]) observableObject() @property
    {
        return _subject;
    }


    Sinks!InputElementType sinks() @property
    {
        typeof(return) dst;
        dst[0] = this;
        return dst;
    }


    Sources!OutputElementType sources() @property
    {
        typeof(return) dst;
        dst[0] = this.observableObject();
        return dst;
    }


    auto subscribe(TObserver)(TObserver observer)
    {
        return doSubscribe(_subject, observer);
    }


  static if(isDuplicable!TImpl)
  {
    RxBlock dup() @property
    {
        RxBlock dst = new RxBlock();
        dst._impl = _impl.dup;
        dst._buffer = _buffer.dup;

        return dst;
    }
  }


  private:
    TImpl _impl;
    SubjectObject!ElementType _subject;
    OutputElementType[] _buffer;


    private this()
    {
        _subject = new SubjectObject!(const(OutputElementType)[])();
    }
}

template makeRxBlock(TImpl)
if(isConverter!TImpl || isHighOrderConverter!TImpl)
{
    auto makeRxBlock(TImpl impl)
    {
        return new RxBlock!TImpl(impl);
    }


    auto makeRxBlock(Args...)(Args args)
    if(is(typeof(TImpl(args))))
    {
        return new RxBlock!TImpl(args);
    }
}


unittest
{
    import rx;
    import std.typecons;

    static struct Impl
    {
        alias InputElementType = int;
        alias OutputElementType = int;

        this(long* ptr) { _sum = ptr; }


        void opCall(in InputElementType[] input, ref OutputElementType[] output)
        {
            import std.algorithm : sum;

            *_sum += input.sum();
        }


        long sum() @property { return *_sum; }
        void sum(long v) @property { *_sum = v;  }

        long* sumPtr() @property { return _sum; }

      private:
        long* _sum;

        long psum() @property { return *_sum; }
    }


    {
        long* ptr = new long;
        auto block = makeRxBlock!Impl(ptr);
        
        auto disp = iota(1024).array.chunks(24).asObservable().doSubscribe(block);
        scope(exit) disp.dispose();

        assert(*ptr == 1023*512);
        assert(block["sum"] == 1023 * 512);
        assert(block["sumPtr"] == ptr);
    }
    {
        long* ptr1 = new long;
        long* ptr2 = new long;
        auto block1 = makeRxBlock!Impl(ptr1);
        auto block2 = makeRxBlock!Impl(ptr2);

        block1.doSubscribe(block2);
    }
}

unittest
{
    import rx;
    import std.typecons;
    import std.algorithm;

    static struct Impl
    {
        alias InputElementType = int;
        alias OutputElementType = int;

        this(long* ptr) { _sum = ptr; }

        void opCall(in InputElementType[] input, ref OutputElementType[] output, scope void delegate() flusher)
        {
            foreach(i; 0 .. 2){
                output.length = input.length;
                output[] = input[];
                *_sum += input.sum;
                flusher();
            }
        }


      private:
        long* _sum;
    }

    static assert(isHighOrderConverter!Impl);

    {
        long* ptr = new long;
        auto block = makeRxBlock!Impl(ptr);
        
        auto disp = iota(1024).array.chunks(24).asObservable().doSubscribe(block);
        scope(exit) disp.dispose();

        assert(*ptr == 1023*512*2);
    }
    {
        long* ptr1 = new long;
        long* ptr2 = new long;
        auto block1 = makeRxBlock!Impl(ptr1);
        auto block2 = makeRxBlock!Impl(ptr2);

        block1.doSubscribe(block2);
    }
}


template RxAdapter(TImpl)
if(isConverter!TImpl || isHighOrderConverter!TImpl)
{
    import rx;

    struct ObserverAdapter(TObserver)
    {
        mixin SimpleObserverImpl!(TObserver, const(TImpl.InputElementType)[]);


        static if(hasCompleted!TObserver || hasFailure!TObserver)
        {
            this(TObserver observer, Disposable disposable, TImpl impl)
            {
                _observer = observer;
                _disposable = disposable;
                _impl = impl;
            }


            this(Args...)(TObserver observer, Disposable disposable, Args args)
            if(is(typeof(TImpl(args))))
            {
                this(observer, disposable, TImpl(args));
            }
        }
        else
        {
            this(TObserver observer, TImpl impl)
            {
                _observer = observer;
                _impl = impl;
            }


            this(Args...)(TObserver observer, Args args)
            if(is(typeof(TImpl(args))))
            {
                this(observer, TImpl(args));
            }
        }

      private:
        void putImpl(const(TImpl.InputElementType)[] input)
        {
            static if(isHighOrderConverter!TImpl)
            {
                _impl(input, _buffer, (){ .put(_observer, _buffer); });
            }
            else
            {
                _impl(input, _buffer);
                .put(_observer, _buffer);
            }
        }


      private:
        TImpl _impl;
        TImpl.OutputElementType[] _buffer;
    }


    struct ObservableAdapter(TObservable, Args...)
    {
        alias ElementType = const(TImpl.OutputElementType)[];

        this(TObservable observable, Args args)
        {
            _observable = observable;
            _args = args;
        }


        auto subscribe(TObserver)(TObserver observer)
        {
            alias ObserverType = ObserverAdapter!(TObserver);
            static if(hasCompleted!TObserver || hasFailure!TObserver)
            {
                auto disposable = new SingleAssignmentDisposable;
                disposable.setDisposable(disposableObject(doSubscribe(
                        _observable,
                        ObserverType(observer, disposable, _args))));
                
                return disposable;
            }
            else
            {
                return doSubscribe(_observable, ObserverType(observer, _args));
            }
        }


      private:
        TObservable _observable;
        Args _args;
    }
}

unittest
{
    import rx;
    import std.typecons;

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
    auto applyImpl(TObservable)(TObservable observable, long* ptr)
    {
        return RxAdapter!(Impl).ObservableAdapter!(TObservable, long*)(observable, ptr);
    }


    long* ptr = new long;

    long* ptr2 = new long;
    auto block = makeRxBlock!Impl(ptr2);
    auto disp = applyImpl(iota(1024).asObservable(), ptr).doSubscribe(block);
    scope(exit) disp.dispose();

    assert(*ptr == 1023*512);
    assert(*ptr2 == 1023*512);
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
        if(!_frontIsComputed){
            cast(void)this.front;
        }

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


  static if(isForwardRange!R && isDuplicable!TImpl)
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
    import rx;
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
        if(!_frontIsComputed){
            cast(void)this.front;
        }

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


  static if(isForwardRange!R && isDuplicable!TImpl)
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
    import rx;
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
        if(!_frontIsComputed){
            cast(void)this.front;
        }

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


  static if(isForwardRange!R && isDuplicable!TImpl)
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
    import rx;
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


template connectTo(TImpl)
{
    import rx : isObservable;

    auto connectTo(S, X...)(S self, X args)
    if(!is(typeof((){ connectTo_old!TImpl(self, args); })))
    {
      static if(isInputRange!S)
        return RangeAdapter!(TImpl, S)(self, args);
      else static if(isObservable!(S, const(TImpl.InputElementType)[]))
        return RxAdapter!(TImpl).ObservableAdapter!(S, X)(self, args);
      else
        static assert(0, S.stringof ~" is not supported by \"connectTo(self, args).\"");
    }
}


auto connectTo(S, TImpl)(S self, TImpl impl)
if(isConverter!TImpl || isHighOrderConverter!TImpl)
{
    return connectTo!TImpl(self, impl);
}


auto connectTo(S, Block)(S self, Block block, ref Disposable disposable)
if(is(typeof(self.doSubscribe(block)) : Disposable))
{
    disposable = self.doSubscribe(block);
    return block;
}


unittest 
{
    import std.range;
    import std.algorithm;

    auto nums = iota(1024);
    auto convd = nums.connectTo!(ArrayConverterOf!(function(int a){ return a + 1; }, int))();
    assert(equal(convd, iota(1, 1025)));
}

unittest 
{
    import std.range;
    import std.algorithm;

    ArrayConverterOf!(function(int a){ return a + 1; }, int) conv;

    auto nums = iota(1024);
    auto convd = nums.connectTo(conv);
    assert(equal(convd, iota(1, 1025)));
}

unittest 
{
    import std.range;
    import std.algorithm;

    auto nums = iota(1024);
    auto convd = nums.connectTo!(ArrayConverterOf!(function(int s, int a){ return a + s; }, int, int))(5);
    assert(equal(convd, iota(5, 1029)));
}

unittest 
{
    import std.range;
    import std.algorithm;

    auto conv = ArrayConverterOf!(function(int s, int a){ return a + s; }, int, int)(5);

    auto nums = iota(1024);
    auto convd = nums.connectTo(conv);
    assert(equal(convd, iota(5, 1029)));
}

unittest 
{
    import std.range;
    import std.algorithm;

    auto conv = ArrayConverterOf!(function(int s, int a){ return a + s; }, int, int)(5);
    
    auto subject = new SubjectObject!int;
    auto block = makeRxBlock(conv);

    int[] results;
    auto nums = iota(1024);
    Disposable[2] disposable;
    disposable[1] = subject.connectTo(block, disposable[0]).doSubscribe!((a){ results ~= a; });
    scope(exit) disposable[0].dispose();
    scope(exit) disposable[1].dispose();

    .put(subject, nums.chunks(64));

    import std.stdio;
    assert(equal(results, iota(5, 1029)));
}


// copy from rx: https://github.com/lempiji/rx/blob/master/source/rx/observer.d
mixin template SimpleObserverImpl(TObserver, E)
{
    import rx;

  public:
    void put(E obj)
    {
        static if (hasFailure!TObserver)
        {
            try
            {
                putImpl(obj);
            }
            catch (Exception e)
            {
                _observer.failure(e);
                _disposable.dispose();
            }
        }
        else
        {
            putImpl(obj);
        }
    }

    static if (hasCompleted!TObserver)
    {
        void completed()
        {
            _observer.completed();
            _disposable.dispose();
        }
    }
    static if (hasFailure!TObserver)
    {
        void failure(Exception e)
        {
            _observer.failure(e);
            _disposable.dispose();
        }
    }
  private:
    TObserver _observer;
    static if (hasCompleted!TObserver || hasFailure!TObserver)
    {
        Disposable _disposable;
    }
}


enum isDuplicable(Conv) = is(typeof((Conv conv){
    Conv other = conv.dup;
}));


/**
1要素を変換するコンバータを配列コンバータへ変換します
*/
struct ArrayConverterOf(Conv)
if(isOneElementConverter!Conv)
{
    alias InputElementType = Conv.InputElementType;
    alias OutputElementType = Conv.OutputElementType;


    this(Conv conv)
    {
        _impl = conv;
    }


    this(X...)(X args)
    if(is(typeof(Conv(args))))
    {
        _impl = Conv(args);
    }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        _impl(input, output);
    }


    void opCall(in InputElementType[] input, ref OutputElementType[] output)
    {
        output.length = input.length;
        foreach(i; 0 .. input.length)
            _impl(input[i], output[i]);
    }


    void opCall(R)(R input, ref OutputElementType[] output)
    if(isInputRange!R && !hasLength!R && is(ElementType!R : InputElementType))
    {
        size_t i;
        while(!input.empty){
            if(output.length < i + 1) output.length = i + 1;
            _impl(input.front, output[i]);
            input.popFront();
        }
    }


  static if(is(typeof((const Conv  conv){ Conv other = conv.dup; })))
  {
    ArrayConverterOf dup() const
    {
        return ArrayConverterOf(_impl.dup);
    }
  }
  else static if(isDuplicable!Conv)
  {
    ArrayConverterOf dup()
    {
        return ArrayConverterOf(_impl.dup);
    }
  }


  private:
    Conv _impl;
}

unittest 
{
    static struct AddConst
    {
        alias InputElementType = int;
        alias OutputElementType = int;

        this(int c) { _c = c; }

        void opCall(InputElementType input, ref OutputElementType output)
        {
            output = input + _c;
        }

      private:
        int _c;
    }

    auto conv = ArrayConverterOf!AddConst(AddConst(2));
    auto input = new int[10];
    int[] output;
    conv(input, output);
    assert(output.length == 10);
    foreach(e; output) assert(e == 2);
}


/**
1要素を変換する関数を配列コンバータへ変換します
*/
struct ArrayConverterOf(alias fn, E, State...)
if(is(typeof(fn(State.init, E.init))))
{
    alias InputElementType = E;
    alias OutputElementType = typeof(fn(State.init, E.init));


  static if(State.length >= 1)
  {
    this(State state)
    {
        _state = state;
    }
  }


    void opCall(InputElementType input, ref OutputElementType output)
    {
        output = fn(_state, input);
    }


    void opCall(R)(R input, ref OutputElementType[] output)
    if(isInputRange!R && hasLength!R && is(ElementType!R : InputElementType))
    {
        output.length = input.length;

        foreach(i; 0 .. input.length)
            output[i] = fn(_state, input[i]);
    }


    void opCall(R)(R input, ref OutputElementType[] output)
    if(isInputRange!R && !hasLength!R && is(ElementType!R : InputElementType))
    {
        size_t i;
        while(!input.empty){
            if(output.length < i + 1) output.length = i + 1;
            output[i] = fn(_state, input.front);
            input.popFront();
        }
    }


  private:
    State _state;
}

unittest 
{
    alias TypeOfInput = int;
    alias State = int;

    auto conv = ArrayConverterOf!((state, input) => input + state, TypeOfInput, State)(2);
    auto input = new int[10];
    int[] output;
    conv(input, output);
    assert(output.length == 10);
    foreach(e; output) assert(e == 2);
}


struct Rechunker(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(size_t size)
    {
        _buffer.length = size;
    }


    void opCall(const(C)[] input, ref C[] output, scope void delegate() flusher)
    {
        immutable len = _buffer.length;

        if(_pos + input.length >= len)
        {
            output.length = len;
            output[0 .. _pos] = _buffer[0 .. _pos];
            output[_pos .. $] = input[0 .. len - _pos];
            flusher();

            input = input[len - _pos .. $];
            while(input.length >= len){
                output[] = input[0 .. len];
                flusher();
                input = input[len .. $];
            }

            _buffer[0 .. input.length] = input[];
            _pos = input.length;
        }
        else
        {
            _buffer[_pos .. _pos + input.length] = input[];
            _pos += input.length;
        }
    }


  private:
    C[] _buffer;
    size_t _pos;
}

unittest 
{
    import std.algorithm : map;
    import rx;
    auto rechunker32 = Rechunker!int(32);

    bool ok = false;
    iota(1024).chain(iota(1024).retro).map!(n => new int[n])
    .asObservable.connectTo(rechunker32).doSubscribe!((a){ assert(a.length == 32); ok = true; });

    assert(ok);
}


struct Dropper(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;


    this(size_t size)
    {
        _rem = size;
    }


    void opCall(const(C)[] input, ref C[] output)
    {
        if(_rem >= input.length){
            _rem -= input.length;
            input.length = 0;
        }else if(_rem > 0){
            input = input[_rem .. $];
            _rem = 0;
        }

        output.length = input.length;
        output[] = input[];
    }


  private:
    size_t _rem;
}

unittest 
{
    import std.range;
    import std.algorithm;

    bool ok = false;
    auto r = iota(2048).array.chunks(32).connectTo(Dropper!int(1024)).joiner;
    assert(equal(r, iota(1024, 2048)));
}


struct NopConverter(C)
{
    alias InputElementType = C;
    alias OutputElementType = C;

    void opCall(InputElementType input, ref OutputElementType output) const
    {
        // output = input + _c;
        output = input;
    }


    void opCall(R)(R input, ref OutputElementType[] output) const
    if(hasLength!R)
    {
        output.length = input.length;

        static if(is(R : const(InputElementType)[])){
            output[] = input[];
        }else{
            foreach(i; 0 .. output.length){
                output[i] = input.front;
                input.popFront();
            }
        }
    }



    NopConverter dup() const @property { return this; }
}

unittest
{
    NopConverter!int c;

    int[] arr;
    c([1, 2, 3], arr);
    assert(arr == [1, 2, 3]);
}


class SerialRxBlocks(I, bool isDuplicable):
    ISPBlock!(I, I, isDuplicable),
    Observer!(const(I)[])
{
    import std.ascii;
    import std.algorithm;
    import std.conv;

    alias InputElementType = I;
    alias OutputElementType = I;
    alias ElementType = I;


    this(ISPBlock!(I, I, isDuplicable)[] blocks...)
    {
        _blocks = blocks.dup;
        foreach(i; 0 .. _blocks.length - 1)
            this._blocks[i].observableObject.doSubscribe(this._blocks[i + 1]);
    }


    bool hasGetter(string key)
    {
        try {
            if(key.front.isDigit){
                auto ns = key.findSplit(".");
                if(ns[1].empty) return false;
                auto n = to!size_t(ns[0]);
                if(n >= _blocks.length) return false;
                return _blocks[n].hasGetter(ns[2]);
            }
        }catch(Exception ex){
            return false;
        }

        return false;
    }


    bool hasSetter(string key)
    {
        try {
            if(key.front.isDigit){
                auto ns = key.findSplit(".");
                if(ns[1].empty) return false;
                auto n = to!size_t(ns[0]);
                if(n >= _blocks.length) return false;
                return _blocks[n].hasSetter(ns[2]);
            }
        }catch(Exception ex){
            return false;
        }

        return false;
    }



    Variant opIndex(string key)
    {
        if(key.front.isDigit){
            auto ns = key.findSplit(".");
            if(!ns[1].empty){
                auto n = to!size_t(ns[0]);
                return _blocks[n][ns[2]];
            }
        }

        {
            import core.exception;
            onRangeError();

            Variant null_;
            return null_;
        }
    }


    void opIndexAssignImpl(Variant value, string key)
    {
        if(key.front.isDigit){
            auto ns = key.findSplit(".");
            if(!ns[1].empty){
                auto n = to!size_t(ns[0]);
                _blocks[n].opIndexAssignImpl(value, ns[2]);
                return;
            }
        }

        {
            import core.exception;
            onRangeError();
        }
    }


    void opIndexAssign(T)(T value, string key)
    {
        static if(is(T : Variant))
            this.opIndexAssignImpl(value, key);
        else
            this.opIndexAssignImpl(Variant(value), key);
    }


    void put(const(I)[] input)
    {
        _blocks[0].put(input);
    }


    void completed() {}
    void failure(Exception ex) {}


    Observable!(const(I)[]) observableObject() @property
    {
        return _blocks[$-1].observableObject;
    }


    Sinks!InputElementType sinks() @property
    {
        return _blocks[0].sinks;
    }


    Sources!OutputElementType sources() @property
    {
        return _blocks[$-1].sources;
    }


    auto subscribe(TObserver)(TObserver observer)
    {
        return doSubscribe(_block[$-1], observer);
    }


  static if(isDuplicable)
  {
    SerialRxBlocks dup() @property
    {
        ISPBlock!(I, I, isDuplicable)[] bs;
        foreach(e; _blocks)
            bs ~= cast(ISPBlock!(I, I, isDuplicable))e.dup;

        return new SerialRxBlocks(bs);
    }
  }


  private:
    ISPBlock!(I, I, isDuplicable)[] _blocks;
}