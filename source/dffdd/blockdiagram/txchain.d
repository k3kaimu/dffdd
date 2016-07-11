module dffdd.blockdiagram.txchain;


struct TXChain
{
    static 
    TXChainImpl!R makeBlock(R)(R r, real alphaIQ, real phiIQ, Gain gain, Voltage iip3, Voltage iip5, Voltage iip7)
    {
        return TXChainImpl!R(r, alphaIQ, phiIQ, gain, iip3, iip5, iip7);
    }


    // 多項係数
    static
    uint mc(uint p, uint k1, uint k2, uint k3)
    {
        uint v = 1;
        foreach(i; 1 .. p+1) v *= i;
        foreach(i; 1 .. k1+1) v /= i;
        foreach(i; 1 .. k2+1) v /= i;
        foreach(i; 1 .. k3+1) v /= i;

        return v;
    }

    unittest
    {
        assert(mc(1, 1, 0, 0) == 1);
        assert(mc(1, 0, 1, 0) == 1);
        assert(mc(1, 0, 0, 1) == 1);

        assert(mc(2, 1, 1, 0) == 2);

        assert(mc(4, 3, 1, 0) == 4);
        assert(mc(4, 2, 2, 0) == 6);
        assert(mc(4, 2, 1, 1) == 12);
    }


    static
    real coefOfABSPower(uint p, real alphaIQ)
    {
        p = p / 2;

        real sum = 0;
        foreach(k2; 0 .. p+1)
            if((p - k2)  % 2 == 0)
                sum += mc(p, (p-k2)/2, k2, (p-k2)/2) * alpha^^(p - k2) * (1 + alphaIQ^^2)^^k2;

        return sum;
    }

    unittest
    {
        assert(coefOfABSPower(100, 0) == 0);
        assert(coefOfABSPower(3-1, 1) == 2);
    }


    static
    struct TXChainImpl(R)
    {
        alias C = Unqual!(ElementType!R);
        alias F = typeof(C.init.re);

        this(R)(R r, real alphaIQ, real phiIQ, Gain gain, Voltage iip3, Voltage iip5, Voltage iip7)
        {
            _r = r;

            _gainQ = alphaIQ * ComplexMethods!C.expi(phiIQ);

            _gain1 = gain.gain;

            _gain3 = (iip3.V == 0 ? 0 : (gain.gain / iip3.V^^2)) * coefOfABSPower(3-1, alphaIQ);
            _gain5 = (iip3.V == 0 ? 0 : (gain.gain / iip3.V^^4)) * coefOfABSPower(5-1, alphaIQ);
            _gain7 = (iip3.V == 0 ? 0 : (gain.gain / iip3.V^^6)) * coefOfABSPower(7-1, alphaIQ);
        }


        auto front() @property
        {
            auto x = _r.front;
            auto x1p = x.re^^2 + x.im^^2,
                 x3 = x * x1p,
                 x5 = x3 * x1p;
                 x7 = x5 * x1p;

            auto xr = x * _gain1 + x3 * _gain3 + x5 * _gain5 + x7 * _gain7;
            return xr + xr.conj * _gainQ;
        }


        void popFront()
        {
            _r.popFront();
        }


        bool empty() @property
        {
            return _r.empty();
        }


      private:
        R _r;
        C _gainQ;
        F _gain1;
        F _gain3;
        F _gain5;
        F _gain7;
    }
}