module app;

import std.algorithm;
import std.range;
import std.datetime;
import std.stdio;
import std.random;
import std.complex;
import std.math;

import dffdd.utils.unit;
import dffdd.blockdiagram.utils;

enum N = 1024*1024;


alias C = Complex!float;

ulong[] rnds;
C[] signal;


void main()
{
    foreach(i; 0 .. N){
        rnds ~= rndGen.front;
        rndGen.popFront();
        signal ~= C(uniform01(), uniform01());
    }

    // writeln("Fiber");
    // benchmarkFiber();
    foreach(i; 0 .. 10){
        realisticBenchmarkHandmade();
    }

    writeln("Handmade");
    realisticBenchmarkHandmade();

    writeln("Buffered Rx");
    // benchmarkRxBuffered();
    realisticBenchmarkRxBuffered();

    writeln("Buffered Range");
    realisticBenchmarkRangeBuffered();

    writeln("Range");
    // benchmarkRange();
    realisticBenchmarkRange();

    writeln("Rx");
    // benchmarkRx();
}


void benchmarkFiber()
{
    int sum;
    auto inst = makeInstrument!(long, (r){ sum = cast(int)r.fold!"a+b"(0L); });
    auto nums = rnds;
    
    StopWatch sw;
    sw.start();
    .put(inst, nums);
    sw.stop();

    writefln("%s", sum);
    writefln("Time: %s[ms]", sw.peek.msecs);
    writefln("Rate: %s[k samples/sec]", N / sw.peek.msecs);
}


void benchmarkRange()
{
    int sum;
    auto nums = rnds;

    StopWatch sw;
    sw.start();
    sum = cast(int)nums.inputRangeObject.fold!"a+b"(0L);
    sw.stop();

    writefln("%s", sum);
    writefln("Time: %s[ms]", sw.peek.msecs);
    writefln("Rate: %s[k samples/sec]", N / (sw.peek.msecs));
}


void benchmarkRx()
{
    import rx;

    int sum;
    auto nums = rnds;

    auto subject = new SubjectObject!long;

    StopWatch sw;
    sw.start();
    auto disposable = subject.doSubscribe!((a){ sum += a; });
    scope(exit) disposable.dispose();
    // .put(subject, nums);

    nums.asObservable().doSubscribe(subject);

    sw.stop();

    writefln("%s", sum);
    writefln("Time: %s[ms]", sw.peek.msecs);
    writefln("Rate: %s[k samples/sec]", N / (sw.peek.msecs));
}


void benchmarkRxBuffered()
{
    import rx;

    int ss;
    auto nums = rnds;
    
    auto subject = new SubjectObject!(ulong[]);

    StopWatch sw;
    sw.start();
    auto disposable = subject.doSubscribe!((a){ ss += a.sum(); });
    scope(exit) disposable.dispose();

    .put(subject, nums.chunks(64));
    sw.stop();

    writefln("%s", ss);
    writefln("Time: %s[ms]", sw.peek.msecs);
    writefln("Rate: %s[k samples/sec]", N / (sw.peek.msecs));
}



void realisticBenchmarkRange()
{
    import dffdd.blockdiagram.iqmixer;
    import dffdd.blockdiagram.amplifier;
    import dffdd.blockdiagram.filter;
    import std.range;

    C[] taps;
    foreach(i; 0 .. 64)
        taps ~= C(0.01, 0.01);

    StopWatch sw;
    sw.start();
    auto ss = signal
        .connectTo(makeIQImbalancer!C(0.dB, (20).dB, 0))
        .connectTo(makeRappModel!C(0.dB, 1, 1))
        .connectTo(makeFIRFilter!C(taps))
        .connectTo(makeRappModel!C(0.dB, 3, 1))
        .connectTo(makeIQImbalancer!C(0.dB, (20).dB, 0)).sum();
    sw.stop();

    writefln("%s", ss);
    writefln("Time: %s[ms]", sw.peek.msecs);
    writefln("Rate: %s[k samples/sec]", N / (sw.peek.msecs));
}


void realisticBenchmarkRangeBuffered()
{
    import dffdd.blockdiagram.iqmixer;
    import dffdd.blockdiagram.amplifier;
    import dffdd.blockdiagram.filter;
    import std.range;

    C[] taps;
    foreach(i; 0 .. 64)
        taps ~= C(0.01, 0.01);

    StopWatch sw;
    sw.start();
    auto ss = signal.chunks(320)
        .connectTo(makeIQImbalancer!C(0.dB, (20).dB, 0))
        .connectTo(makeRappModel!C(0.dB, 1, 1))
        .connectTo(makeFIRFilter!C(taps))
        .connectTo(makeRappModel!C(0.dB, 3, 1))
        .connectTo(makeIQImbalancer!C(0.dB, (20).dB, 0)).map!(a => a.sum()).sum();
    sw.stop();

    writefln("%s", ss);
    writefln("Time: %s[ms]", sw.peek.msecs);
    writefln("Rate: %s[k samples/sec]", N / (sw.peek.msecs));
}


void realisticBenchmarkRxBuffered()
{
    import rx;

    import dffdd.blockdiagram.iqmixer;
    import dffdd.blockdiagram.amplifier;
    import dffdd.blockdiagram.filter;
    import std.range;

    C[] taps;
    foreach(i; 0 .. 64)
        taps ~= C(0.01, 0.01);

    C ss = C(0, 0);

    StopWatch sw;
    sw.start();
    signal.chunks(320).asObservable
        .connectTo(makeIQImbalancer!C(0.dB, (20).dB, 0))
        .connectTo(makeRappModel!C(0.dB, 1, 1))
        .connectTo(makeFIRFilter!C(taps))
        .connectTo(makeRappModel!C(0.dB, 3, 1))
        .connectTo(makeIQImbalancer!C(0.dB, (20).dB, 0))
        .doSubscribe!((a){ ss += a.sum(); });
    sw.stop();

    writefln("%s", ss);
    writefln("Time: %s[ms]", sw.peek.msecs);
    writefln("Rate: %s[k samples/sec]", N / (sw.peek.msecs));
}



void realisticBenchmarkHandmade()
{
    import dffdd.blockdiagram.filter;

    C[] taps;
    foreach(i; 0 .. 64)
        taps ~= C(0.01, 0.01);


    StopWatch sw;
    sw.start();
    auto input = signal.dup;

    foreach(ref e; input)
        e = e + e.conj * 0.1;

    
    foreach(ref e; input) {
        auto x = e;
        auto r = sqAbs(x);

        r = 1.0 / (sqrt( 1 + r));
        e = r * x;
    }

    // input = input.chunks(1024).connectTo(makeFIRFilter(taps)).joiner.map!(a => C(a.re, a.im)).array();
    {
        auto recv = input.chunks(320).connectTo(makeFIRFilter(taps));
        auto dst = input;
        while(!recv.empty){
            auto front = recv.front;
            recv.popFront();

            dst[0 .. front.length] = front[];
            dst = dst[front.length .. $];
        }
    }

    foreach(ref e; input) {
        auto x = e;
        auto r = sqAbs(x);

        e = x / (sqrt(1 + r^^3).cbrt);
    }

    foreach(ref e; input)
        e = e + e.conj * 0.1;

    C ss = input.sum();

    sw.stop();

    writefln("%s", ss);
    writefln("Time: %s[ms]", sw.peek.msecs);
    writefln("Rate: %s[k samples/sec]", N / (sw.peek.msecs));
}


unittest
{
    alias C = Complex!float;

    auto source = new SubjectObject!C;
    auto iqi = makeRxBlock(makeIQImbalancer!C(0.dB, (20).dB, 0));
    auto rapp = makeRxBlock(makeRappModel!C(0.dB, 1, 1));
    auto sink = new SubjectObject!C;

    source.doSubscribe(iqi);
    iqi.doSubscribe(rapp);
    rapp.doSubscribe(sink);
}

unittest
{
    alias C = Complex!float;

    C[] signal;
    auto tx = signal.chunks(320)
        .connectTo(makeIQImbalancer!C(0.dB, (20).dB, 0))
        .conenctTo(makeRappModel!C(0.dB, 1, 1));
}

unittest
{
    alias C = Complex!float;

    auto subject = new SubjectObject!C;
    auto tx = subject
        .connectTo(makeIQImbalancer!C(0.dB, (20).dB, 0))
        .connectTo(makeRappModel!C(0.dB, 1, 1));

}