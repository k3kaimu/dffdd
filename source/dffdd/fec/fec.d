module dffdd.fec.fec;

import std.math;
import std.complex;

import dffdd.mod.primitives : Bit;

interface Fec(T)
{
    Bit[] encode(in Bit[] info, return ref Bit[] encoded);
    Bit[] decode(in T[] input, return ref Bit[] decoded);
}


interface BlockFec(T) : Fec!(T)
{
    size_t inputLength() const;
    size_t outputLength() const;


    override
    Bit[] encode(in Bit[] info, return ref Bit[] encoded)
    in {
        assert(info.length % this.inputLength == 0);
    }
    out(r) {
        // assert(r.length % this.outputLength == 0);
        // assert(info.length / this.inputLength == r.length / this.outputLength);
    };


    override
    Bit[] decode(in T[] input, return ref Bit[] decoded)
    in {
        assert(input.length % this.outputLength == 0);
    }
    out(r) {
        // assert(r.length % this.inputLength == 0);
        // assert(input.length / this.outputLength == r.length / this.inputLength);
    };
}


struct P0P1Calculator(C, F)
{
    this(Mod)(Mod mod)
    in(mod.symOutputLength == 1)
    {
        _nbits = mod.symInputLength;
        _p0 = new F[_nbits];
        _p1 = new F[_nbits];

        C[] sym = new C[1];
        Bit[] bits = new Bit[_nbits];
        foreach(n; 0 .. 1 << _nbits) {
            foreach(k; 0 .. _nbits)
                bits[k] = (n & (1 << k)) ? 1 : 0;

            mod.modulate(bits, sym);
            _points ~= sym[0];
        }
    }


    ref F[] computeP0P1(in C[] inputs, return ref F[] p0p1, F N0)
    {
        p0p1.length = inputs.length * _nbits;

        foreach(i, x; inputs) {
            _p0[] = 0;
            _p1[] = 0;

            foreach(j, y; _points) {
                auto p = exp(-(x - y).sqAbs/ N0);

                foreach(k; 0 .. _nbits) {
                    if(j & (1 << k))
                        _p1[k] += p;
                    else
                        _p0[k] += p;
                }
            }

            foreach(k; 0 .. _nbits)
                p0p1[i * _nbits + k] = _p0[k] / _p1[k];
        }

        return p0p1;
    }


  private:
    size_t _nbits;
    C[] _points;
    F[] _p0;
    F[] _p1;
}


auto p0p1Calculator(F, Mod)(Mod mod)
{
    return P0P1Calculator!(Mod.OutputElementType, F)(mod);
}
