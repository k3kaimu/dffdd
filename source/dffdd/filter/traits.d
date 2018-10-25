module dffdd.filter.traits;

import std.traits;

enum bool isBlockConverter(T, A, B) = is(typeof((T converter){
    A[] input;
    B[] output;
    converter(input, output);
    size_t lenInput = converter.inputBlockLength;
    // size_t dimOutput = converter.outputDim;
}));


/**
+ opCallImpl(A, ref B)
が定義されている時，次の定義を追加する
*/
mixin template ConverterOpCalls(A, B)
{
    void opCall(A a, ref B b)
    {
        this.opCallImpl(a, b);
    }


    B opCall(A a)
    {
        typeof(return) b;
        this.opCallImpl(a, b);
        return b;
    }


    void opCall(A[] as, ref B[] bs)
    {
        if(bs.length != as.length) bs.length = as.length;
        foreach(i; 0 .. as.length)
            this.opCallImpl(as[i], bs[i]);
    }


    B[] opCall(A)(A[] as)
    {
        typeof(return) bs;
        this.opCall(as, bs);
        return bs;
    }
}


enum bool isFilter(F, A, B) = is(typeof((F filter){
    A[] input;
    B[] desires;
    B[] errors;
    size_t lenInput = filter.inputBlockLength;
    filter.apply!(Yes.learning)(input, desires, errors);
    filter.apply!(No.learning)(input, desires, errors);
}));


enum bool isState(S) = is(typeof((ref S s){
    auto ssp = &(s.state);
    static assert(isStaticArray!(typeof(*ssp)));

    auto sep = &(ssp[0]);
    static assert(isStaticArray!(typeof(*sep)));

    alias C = typeof(s.state[0][0]);

    auto wsp = &(s.weight);
    static assert(isStaticArray!(typeof(*wsp)));
    static assert(typeof(*wsp).length == typeof(*ssp).length);

    auto wep = &(wsp[0]);
    static assert(isStaticArray!(typeof(*wep)));
    static assert(typeof(*wep).length == typeof(*sep).length);

    C c;
    s.update(c);

    // .power is option
    //auto psp = &(s.power);
    //static assert(isStaticArray!(typeof(*psp)));
    //static assert(typeof(*psp).length == typeof(*sep).length);

    C v = s.value;
}));


enum bool isWeightAdapter(U, S) = is(typeof((ref U u, ref S s){
    enum bool usePower = U.usePower;

    alias C = typeof(s.state[0][0]);
    C e;

    u.adapt(s, e);
}));


unittest
{
    import std.complex;

    alias C = Complex!float;

    static struct TestState1
    {
        C[1][1] weight;
        C[1][1] state;
        float[1] power;
        C _value;

        void update(C c)
        {
            state[0][0] = c; power[0] = c.abs()^^2;
        }

        C value() { return weight[0][0] * state[0][0]; }
    }

    static assert(isState!(TestState1));

    static struct TestAdapter1
    {
        enum bool usePower = false;

        void adapt(ref TestState1 s, C e)
        {
            s.weight[0][0] += e;
        }
    }

    static assert(isWeightAdapter!(TestAdapter1, TestState1));
}


enum bool isAdaptiveFilter(F, C = Complex!float) = is(typeof((F f, C[] buf){
    f.apply!false(cast(const)buf, cast(const)buf, buf);
    f.apply!true(cast(const)buf, cast(const)buf, buf);
}));
