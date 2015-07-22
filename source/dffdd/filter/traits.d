module dffdd.filter.traits;

import std.traits;

enum bool isState(S) = is(typeof((ref S s){
    auto ssp = &(s.state);
    static assert(isStaticArray!(typeof(*ssp)));

    auto sep = &(ss[0]);
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

        void update(C c)
        {
            state[0][0] = c; power[0] = c.abs()^^2;
        }
    }

    static assert(isState!(TestState1));

    static struct TestAdapter1
    {
        void adapt(ref TestState1 s, C e)
        {
            s.weight[0] += e;
        }
    }

    static assert(isWeightAdapter!(TestAdapter1, TestState1));
}

