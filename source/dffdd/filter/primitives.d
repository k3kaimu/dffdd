module dffdd.filter.primitives;

import std.range;
import std.algorithm;
import std.traits;

import dffdd.filter.traits;
import dffdd.mod.primitives;
import dffdd.utils.fft;


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
