module dffdd.filter.primitives;

import std.range;
import std.algorithm;
import std.traits;
import std.complex;

import dffdd.filter.traits;
import dffdd.mod.primitives;
import dffdd.utils.fft;
import dffdd.filter.state;
import dffdd.math.math;


final class LimitedTrainingAdaptor(Adaptor)
{
    Adaptor parent;
    alias parent this;

    this(Adaptor adaptor, size_t limitSamples)
    {
        parent = adaptor;
        _updateLimite = limitSamples;
        _remain = limitSamples;
    }


    void adapt(State, C)(ref State state, C error)
    {
        if(_remain != 0) {
            parent.adapt(state, error);
            --_remain;
        }
    }


    void resetLimit()
    {
        _remain = _updateLimite;
    }


  private:
    immutable size_t _updateLimite;
    size_t _remain;
}


LimitedTrainingAdaptor!Adaptor trainingLimit(Adaptor)(Adaptor adaptor, size_t limitSamples)
{
    return new typeof(return)(adaptor, limitSamples);
}


final class IgnoreHeadSamplesAdaptor(Adaptor)
{
    Adaptor parent;
    alias parent this;

    this(Adaptor adaptor, size_t ignoreSamples)
    {
        parent = adaptor;
        _ignoreSamples = ignoreSamples;
        _remain = ignoreSamples;
    }


    void adapt(State, C)(ref State state, C error)
    {
        if(_remain != 0) {
            --_remain;
            return;
        } else {
            parent.adapt(state, error);
        }
    }


    void resetIgnore()
    {
        _remain = _ignoreSamples;
    }


  private:
    immutable size_t _ignoreSamples;
    size_t _remain;
}


IgnoreHeadSamplesAdaptor!Adaptor ignoreHeadSamples(Adaptor)(Adaptor adaptor, size_t ignoreSamples)
{
    return new typeof(return)(adaptor, ignoreSamples);
}


interface IAdapter(C)
{
    void adapt(ref MultiFIRState!C state, C error);
}


final class Distorter(C, funcs...)
{
    enum size_t outputDim = funcs.length;
    enum size_t inputBlockLength = 1;


    this() {}


    void opCallImpl(C input, ref C[] output)
    {
        output.length = outputDim;
        foreach(p, f; funcs)
            output[p] = f(input);
    }

    mixin ConverterOpCalls!(const(C), C[]);
}

unittest
{
    import std.complex;
    import dffdd.filter.traits;

    static assert(isBlockConverter!(Distorter!(Complex!float, x => x, x => x*2), Complex!float, Complex!float[]));
}


final class MemoryDistorter(C, size_t M, funcs...)
{
    enum size_t outputDim = funcs.length;
    enum size_t inputBlockLength = 1;

    this(C[M] init) {
        _state = init;
    }

    void opCallImpl(C input, ref C[] outputs)
    {
        outputs.length = outputDim;
        foreach_reverse(j; 1 .. _state.length) _state[j] = _state[j-1];
        _state[0] = input;

        foreach(p, f; funcs)
            outputs[p] = f(_state);
    }

    mixin ConverterOpCalls!(const(C), C[]);

  private:
    C[M] _state;
}


unittest
{
    auto dist1 = new MemoryDistorter!(int, 2, (ref s) => s[0] * s[1])([0, 0]);
    assert(dist1(1)[0] == 0);        // state = [1, 0]
    assert(dist1(2)[0] == 2);        // state = [2, 1]
    assert(dist1(3)[0] == 6);        // state = [3, 2]

    auto dist2 = new MemoryDistorter!(int, 2, (ref s) => s[0])([0, 0]);
    assert(dist2(1)[0] == 1);        // state = [1, 0]
    assert(dist2(2)[0] == 2);        // state = [2, 1]
    assert(dist2(3)[0] == 3);        // state = [3, 2]
}


final class OnlyPADistorter(C, size_t P)
if(P % 2 == 1)
{
    enum size_t outputDim = (P + 1)/2;
    enum size_t inputBlockLength = 1;

    void opCallImpl(C input, ref C[] output)
    {
        output.length = outputDim;
        foreach(p; 1 .. P+1){
            if(p % 2 == 1)
                output[(p-1)/2] = input * (input.sqAbs() ^^ ((p-1)/2));
        }
    }

    mixin ConverterOpCalls!(const(C), C[]);
}


final class SubSetOfPADistorter(C, size_t P)
if(P % 2 == 1)
{
    enum size_t outputDim = P + 1;
    enum size_t inputBlockLength = 1;


    ptrdiff_t indexOfConjugated(size_t i)
    {
        if(i % 2 == 0)
            return i + 1;
        else
            return i - 1;
    }


    this() {}


    void opCallImpl(C input, ref C[] output)
    {
        output.length = outputDim;
        foreach(p; 0 .. P+1){
            if(p % 2 == 0){
                output[p] = input * (input.sqAbs() ^^ (p/2));
            }else{
                output[p] = input.conj * (input.sqAbs() ^^ (p/2));
            }
        }
    }

    mixin ConverterOpCalls!(const(C), C[]);
}

unittest
{
    import std.math;
    import dffdd.filter.traits;
    import std.stdio;

    static assert(isBlockConverter!(SubSetOfPADistorter!(Complex!float, 3), Complex!float, Complex!float[]));

    auto dist1 = new SubSetOfPADistorter!(Complex!float, 1);
    Complex!float[] outputs;
    dist1(Complex!float(1, 1), outputs);
    foreach(i, e; [Complex!float(1, 1), Complex!float(1, -1)])
    {
        assert(isClose(e.re, outputs[i].re));
        assert(isClose(e.im, outputs[i].im));
    }

    auto dist3 = new SubSetOfPADistorter!(Complex!float, 3);
    dist3(Complex!float(1, 1), outputs);
    foreach(i, e; [Complex!float(1, 1), Complex!float(1, -1), Complex!float(2, 2), Complex!float(2, -2)])
    {
        assert(isClose(e.re, outputs[i].re));
        assert(isClose(e.im, outputs[i].im));
    }

    auto dist5 = new SubSetOfPADistorter!(Complex!float, 5);
    dist5(Complex!float(1, 1), outputs);
    foreach(i, e; [Complex!float(1, 1), Complex!float(1, -1), Complex!float(2, 2), Complex!float(2, -2), Complex!float(4, 4), Complex!float(4, -4)])
    {
        assert(isClose(e.re, outputs[i].re));
        assert(isClose(e.im, outputs[i].im));
    }
}


final class PADistorter(C, size_t P)
if(P % 2 == 1)
{
    enum size_t outputDim = (P+1) * ((P+1)/2 + 1) / 2;
    enum size_t inputBlockLength = 1;


    // template indexOfPQ(size_t i, size_t pq = 1)
    // {
    //     enum size_t odim = (pq+1) * ((pq+1)/2 + 1) / 2
    //     static if(i < odim)
    //         enum ptrdiff_t[2] indexOfPQ = indexOfPQ_FindPair!(odim, i - (pq-1) * ((p1-1)/2 + 1) / 2);
    //     else
    //         enum ptrdiff_t[2] indexOfPQ = indexOfPQ!(i, pq+2);
    // }


    // template indexOfPQ_FindPair(size_t pq, size_t j)
    // {
    //     enum ptrdiff_t[2] indexOfPQ_FindPair = [pq-j, j];
    // }


    // template indexOfConjugated(size_t i)
    // {
    //     enum ptrdiff_t[2] pq = indexOfPQ(i);
    //     enum ptrdiff_t M = (pq[0] + pq[1] - 2);

    //     alias indexOfConjugated = (M+1)*((M+1)/2 + 1)/2 + pq[1];
    // }

    static @nogc pure nothrow @safe
    ptrdiff_t indexOfConjugated(size_t k)
    {
        foreach(ptrdiff_t p; 1 .. P+1){
            if(p % 2 == 0) continue;

            auto odim = (p+1) * ((p+1)/2 + 1) / 2;
            if(k < odim) {
                ptrdiff_t odim1 = (p-1) * ((p-1)/2 + 1) / 2;
                ptrdiff_t j = k - odim1;
                ptrdiff_t i = p - j;
                return odim1 + i;
            }
        }

        return -1;
    }


    this() {}


    void opCallImpl(C input, ref C[] output)
    {
        output.length = outputDim;
        size_t i = 0;
        foreach(p; 1 .. P+1){
            if(p % 2 == 0) continue;

            foreach(j; 0 .. p+1){
                output[i] = input^^(p-j) * input.conj^^(j);
                ++i;
            }
        }
    }

    mixin ConverterOpCalls!(const(C), C[]);
}

unittest
{
    import std.math;
    import dffdd.filter.traits;
    import std.stdio;


    enum expi = (real x){ return cast(Complex!float)std.complex.expi(x); };

    static assert(isBlockConverter!(PADistorter!(Complex!float, 3), Complex!float, Complex!float[]));

    auto dist1 = new PADistorter!(Complex!float, 1);
    Complex!float[] outputs;
    dist1(expi(PI_4), outputs);
    foreach(i, e; [expi(PI_4), expi(-PI_4)])
    {
        assert(isClose(e.re, outputs[i].re));
        assert(isClose(e.im, outputs[i].im));
    }

    auto dist3 = new PADistorter!(Complex!float, 3);
    dist3(expi(PI_4)*2, outputs);
    foreach(i, e; [expi(PI_4)*2, expi(-PI_4)*2,
                   expi(PI_4*3)*8, expi(PI_4)*8, expi(-PI_4)*8, expi(-PI_4*3)*8])
    {
        assert(isClose(e.re, outputs[i].re));
        assert(isClose(e.im, outputs[i].im));
    }

    auto dist5 = new PADistorter!(Complex!float, 5);
    dist5(expi(PI_4)*2, outputs);
    foreach(i, e; [expi(PI_4)*2, expi(-PI_4)*2,
                   expi(PI_4*3)*8, expi(PI_4)*8, expi(-PI_4)*8, expi(-PI_4*3)*8,
                   expi(PI_4*5)*32, expi(PI_4*3)*32, expi(PI_4)*32, expi(-PI_4)*32, expi(-PI_4*3)*32, expi(-PI_4*5)*32])
    {
        assert(isClose(e.re, outputs[i].re));
        assert(isClose(e.im, outputs[i].im));
    }
}
unittest
{
    alias D = PADistorter!(Complex!float, 9);
    assert(D.indexOfConjugated(0) == 1);
    assert(D.indexOfConjugated(1) == 0);
    assert(D.indexOfConjugated(2) == 5);
    assert(D.indexOfConjugated(3) == 4);
    assert(D.indexOfConjugated(4) == 3);
    assert(D.indexOfConjugated(5) == 2);
    assert(D.indexOfConjugated(6) == 11);
    assert(D.indexOfConjugated(7) == 10);
    assert(D.indexOfConjugated(8) == 9);
    assert(D.indexOfConjugated(9) == 8);
    assert(D.indexOfConjugated(10) ==7);
    assert(D.indexOfConjugated(11) == 6);
}


final class PADistorterWithEvenOrder(C, size_t P)
{
    enum size_t outputDim = P + P*(P+1)/2 + 1;
    enum size_t inputBlockLength = 1;


    this() {}


    void opCallImpl(C input, ref C[] output)
    {
        output.length = outputDim;
        size_t i = 0;
        foreach(p; 0 .. P+1){
            foreach(j; 0 .. p+1){
                output[i] = input^^(p-j) * input.conj^^(j);
                ++i;
            }
        }
    }

    mixin ConverterOpCalls!(const(C), C[]);

}


final class Laguerre2DDistorter(C, size_t P)
{
    enum size_t outputDim = P + P*(P+1)/2 + 1;
    enum size_t inputBlockLength = 1;


    this() {}


    void opCallImpl(C input, ref C[] output)
    {
        import dffdd.math.math;

        output.length = outputDim;
        size_t i = 0;
        foreach(int p; 0 .. P+1){
            foreach(int j; 0 .. p+1){
                output[i] = laguerreOBF(input, p-j, j);
                ++i;
            }
        }
    }

    mixin ConverterOpCalls!(const(C), C[]);
}


final class LaguerreOnlyPADistorter(C, size_t P)
{
    enum size_t outputDim = (P + 1)/2;
    enum size_t inputBlockLength = 1;

    this() {}

    void opCallImpl(C input, ref C[] output)
    {
        output.length = outputDim;
        foreach(int i; 0 .. cast(int)(P+1)/2){
            output[i] = laguerreOBF(input, i+1, i);
        }
    }

    mixin ConverterOpCalls!(const(C), C[]);
}


final class ZipDistorter(C, DistList...)
{
    import std.meta : staticMap;
    import std.numeric : lcm;


    this(DistList dists)
    {
        _dists = dists;

        {
            size_t sum;
            foreach(ref e; _dists) sum += e.outputDim;
            _dim = sum;
        }

        {
            size_t dst = _dists[0].inputBlockLength;
            static foreach(i; 1 .. DistList.length) dst = lcm(dst, _dists[i].inputBlockLength);
            _blen = dst;
        }
    }


    size_t outputDim() const
    {
        return _dim;
    }


    size_t inputBlockLength() const
    {
        return _blen;
    }


    void opCallImpl(C input, ref C[] output)
    {
        output.length = this.outputDim;
        size_t cnt = 0;
        foreach(k, E; DistList) {
            immutable dk = _dists[k].outputDim;
            auto buf = output[cnt .. cnt + dk];
            _dists[k](input, buf);
            cnt += dk;
        }
    }


    mixin ConverterOpCalls!(const(C), C[]);

  private:
    DistList _dists;
    size_t _dim, _blen;
}


auto zipDistorter(C, DistList...)(DistList distorters)
{
    return new ZipDistorter!(C, DistList)(distorters);
}


unittest
{
    import std.complex;

    auto dist1 = new Distorter!(Complex!float, x => x, x => x*2)();
    auto dist2 = new Distorter!(Complex!float, x => x)();
    auto zipdist = zipDistorter!(Complex!float)(dist1, dist2);
    static assert(isBlockConverter!(typeof(zipdist), Complex!float, Complex!float[]));

    auto x = Complex!float(2, 0);
    Complex!float[] buf;
    zipdist(x, buf);
    assert(buf.length == 3);
    assert(buf[0] == x);
    assert(buf[1] == x*2);
    assert(buf[2] == x);
}


__EOF__

/**
自己干渉除去フィルタの入出力の静的及び動的情報を管理します．
*/
struct InterfaceInfo
{
    /**
    入出力が周波数軸であることを表します．

    + removeCP: サイクリックプレフィックスが除去されているかどうか
    + numOfFFT: FFTサイズ
    */
    struct Frequency(Flag!"removeCP" _removeCP)
    {
        alias removeCP = _removeCP;
        uint numOfFFT;
    }


    /**
    入出力が時間軸であることを表します．

    + removeCP: サイクリックプレフィックスが除去されているかどうか
    + packed: パッキング(チャンキング)されているかどうか
    + packSize: パッキングされているとき，パック一つのサイズ
    */
    struct Time(Flag!"removeCP" _removeCP, Flag!"packed" _packed)
    {
        alias removeCP = _removeCP;
        alias packed = _packed;

      static if(packed)
      {
        size_t packSize;
      }
    }
}


enum bool isAdaptiveFilter(F, C) = is(typeof((F f){
    C[] cpx;
    f.apply!false(cpx, cpx, cpx);
    f.apply!true(cpx, cpx, cpx);

    auto intype = f.inputInterfaceInfo;
    auto outype = f.outputInterfaceInfo;
}));
