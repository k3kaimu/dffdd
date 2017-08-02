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

enum N = 2_000_000L;


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

    writeln("Range");
    // benchmarkRange();
    realisticBenchmarkRange();

    writeln("Buffered Range");
    realisticBenchmarkRangeBuffered();

    writeln("Rx");
    // benchmarkRx();

    writeln("Buffered Rx");
    // benchmarkRxBuffered();
    realisticBenchmarkRxBuffered();

    writeln("Handmade");
    realisticBenchmarkHandmade();
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

    auto txiq = signal.connectTo!(IQImbalanceConverter!C)(0.dB, (20).dB, 0).inputRangeObject;
    auto txamp = txiq.connectTo!(RappModelConverter!C)(0.dB, 1, 1).inputRangeObject;
    auto recv = txamp.connectTo!(FIRFilterConverter!C)(taps).inputRangeObject;
    auto lna = recv.connectTo!(RappModelConverter!C)(0.dB, 1, 1).inputRangeObject;
    auto rxiq = lna.connectTo!(IQImbalanceConverter!C)(0.dB, (20).dB, 0).inputRangeObject;

    StopWatch sw;
    sw.start();
    auto ss = rxiq.sum;
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

    auto txiq = signal.chunks(320).connectTo!(IQImbalanceConverter!C)(0.dB, (20).dB, 0).inputRangeObject;
    auto txamp = txiq.connectTo!(RappModelConverter!C)(0.dB, 1, 1).inputRangeObject;
    auto recv = txamp.connectTo!(FIRFilterConverter!C)(taps).inputRangeObject;
    auto lna = recv.connectTo!(RappModelConverter!C)(0.dB, 1, 1).inputRangeObject;
    auto rxiq = lna.connectTo!(IQImbalanceConverter!C)(0.dB, (20).dB, 0).inputRangeObject.joiner;

    StopWatch sw;
    sw.start();
    auto ss = rxiq.sum;
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

    auto txiq = RxBlock!(IQImbalanceConverter!C)(0.dB, (20).dB, 0);
    auto txamp = RxBlock!(RappModelConverter!C)(0.dB, 1, 1);
    auto recv = RxBlock!(FIRFilterConverter!C)(taps);
    auto lna = RxBlock!(RappModelConverter!C)(0.dB, 1, 1);
    auto rxiq = RxBlock!(IQImbalanceConverter!C)(0.dB, (20).dB, 0);

    // pragma(msg, typeof(sig).stringof);

    C ss = C(0, 0);
    txiq.doSubscribe(txamp);
    txamp.doSubscribe(recv);
    recv.doSubscribe(lna);
    lna.doSubscribe(rxiq);
    rxiq.doSubscribe!((a){ ss += a.sum; });

    StopWatch sw;
    sw.start();
    .put(txiq, signal.chunks(320));
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


    auto input = signal.dup;

    StopWatch sw;
    sw.start();
    foreach(ref e; input)
        e = e + e.conj * 0.1;

    
    foreach(ref e; input) {
        auto x = e;
        auto r = abs(x),
             u = x / r;     // unit vector

        // rが小さすぎるときに，単位ベクトルが発散するのを防ぐ
        if(r <= 1E-6){
            e = x;
        }
        else{
            r = (r) / (sqrt( 1 + (r)^^2 ));
            e = r * u;
        }
    }

    input = input.chunks(1024).connectTo!(FIRFilterConverter!C)(taps).joiner.array().dup;

    foreach(ref e; input) {
        auto x = e;
        auto r = abs(x),
             u = x / r;     // unit vector

        // rが小さすぎるときに，単位ベクトルが発散するのを防ぐ
        if(r <= 1E-6){
            e = x;
        }
        else{
            r = (r) / (sqrt( 1 + (r)^^2 ));
            e = r * u;
        }
    }

    foreach(ref e; input)
        e = e + e.conj * 0.1;

    C ss = C(0, 0);
    foreach(e; input)
        ss += e;

    sw.stop();

    writefln("%s", ss);
    writefln("Time: %s[ms]", sw.peek.msecs);
    writefln("Rate: %s[k samples/sec]", N / (sw.peek.msecs));
}