module dffdd.utils.binary;

import core.bitop;

import std.algorithm;
import std.random;
import std.range;
import std.traits;

import carbon.stream;


struct Bit
{
    ubyte value;
    alias value this;


    this(long v) pure nothrow @safe @nogc
    {
        this.opAssign(v);
    }


    void opAssign(long v) pure nothrow @safe @nogc
    {
        if(v)
            value = 1;
        else
            value = 0;
    }
}


/**
渡されたレンジの要素のうち、下部nDigitsビットだけをレンジにして返します。
*/
auto toBinaryDigits(bool fromLSB = true, R)(R rng, size_t nDigits)
if(isInputRange!R && isUnsigned!(Unqual!(ElementType!R)))
in{
    assert(nDigits <= (ElementType!R).sizeof * 8);
}
body{
    alias E = Unqual!(ElementType!R);

    static struct Result()
    {
        Bit front() @property pure nothrow @safe const
        {
            return Bit((_front & _mask) != 0 ? 1 : 0);
        }


        bool empty() @property pure nothrow @safe const
        {
            return _mask == 0;
        }


        void popFront()
        {
          static if(fromLSB)
          {
            immutable _maskmask = (1 << (_nDigs - 1)) - 1;
            _mask = ((_mask & _maskmask) << 1) | (_mask >> (_nDigs - 1));

            if(_mask == 1){
                _range.popFront();
                if(!_range.empty)
                    _front = _range.front;
                else
                    _mask = 0;
            }
          }
          else
          {
            immutable _maskmask1 = (1 << _nDigs) - 2;
            immutable _maskmask2 = (1 << _nDigs) - 1;
            _mask = ((_mask & _maskmask1) >> 1) | ((_mask << (_nDigs - 1)) & _maskmask2);

            if(_mask == (1 << (_nDigs - 1))){
                _range.popFront();
                if(!_range.empty)
                    _front = _range.front;
                else
                    _mask = 0;
            }
          }
        }


      static if(isForwardRange!R)
        typeof(this) save() @property
        {
            typeof(this) dst = this;
            dst._range = dst._range.save;
            return dst;
        }


      static if(hasLength!R)
        size_t length() @property
        {
          static if(fromLSB)
          {
            return _range.length * _nDigs - bsr(_mask);
          }
          else
          {
            return _range.length * _nDigs - (_nDigs - 1 - bsr(_mask));
          }
        }


      private:
        R _range;
        E _front;
        E _mask;
        size_t _nDigs;
    }

    
    if(!rng.empty && nDigits){
      static if(fromLSB)
        return Result!()(rng, rng.front, 1, nDigits);
      else
        return Result!()(rng, rng.front, 1 << (nDigits - 1), nDigits);
    }else{
        Result!() dst;
        dst._mask = 0;
        
        return dst;
    }
}
///
unittest{
    uint[] input = [0, 1, 2, 3, 4, 5, 6, 7];

    auto binDigs = toBinaryDigits(input, 3);
    assert(equal(binDigs, [0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1]));

    binDigs = toBinaryDigits(input, 2);
    assert(equal(binDigs, [0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1]));

    binDigs = toBinaryDigits(input, 1);
    assert(equal(binDigs, [0, 1, 0, 1, 0, 1, 0, 1]));

    binDigs = toBinaryDigits(input, 0);
    assert(binDigs.empty);

    binDigs = toBinaryDigits(input[0 .. 3], 32);
    assert(equal(binDigs, chain(0.repeat(32), 1.only, 0.repeat(31), [0, 1], 0.repeat(30))));
}
///
unittest{
    uint[] input = [0, 1, 2, 3, 4, 5, 6, 7];

    auto binDigs = toBinaryDigits!false(input, 3);
    assert(equal(binDigs, [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1]));

    binDigs = toBinaryDigits!false(input, 2);
    assert(equal(binDigs, [0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1]));

    binDigs = toBinaryDigits!false(input, 1);
    assert(equal(binDigs, [0, 1, 0, 1, 0, 1, 0, 1]));

    binDigs = toBinaryDigits!false(input, 0);
    assert(binDigs.empty);
}


/**
渡された疑似乱数生成器を用いてランダムなビット列を生成します
*/
auto randomBits(E = ubyte, Rnd)(Rnd rnd = std.random.rndGen)
if(isUniformRNG!Rnd && isUnsigned!(ElementType!Rnd))
{
    alias E = ElementType!Rnd;
    return rnd.toBinaryDigits(E.sizeof * 8).map!"cast(ubyte)(a ? 1 : 0)";
}

unittest
{
    import std.random;
    auto rndbits = randomBits(Random(0));

    alias E = ElementType!Random;
    auto rnd = Random(0);
    foreach(i; 0 .. 1000) {
        E x = rnd.front;
        foreach(j; 0 .. E.sizeof * 8){
            assert(((x >> j) & 1) == rndbits.front);
            rndbits.popFront();
        }
        rnd.popFront();
    }
}

unittest
{
    auto rndbits = randomBits();
}



void oct2bin(Int)(string oct, Int[] bin, int skiplast, int flip)
in{
    assert(bin.length <= oct.length * 3);
}
body{
    auto tmp = new Int[oct.length * 3];
    {
        auto sink = tmp;
        sink.put(oct.map!"a - '0'"().toBinaryDigits!false(3).map!(a => cast(Int)(a ? -1 : 1))());
    }

    if(skiplast)
        bin[] = tmp[0 .. bin.length];
    else
        bin[] = tmp[$ - bin.length .. $];

    if(flip)
        bin.reverse();
}
unittest{
    auto result = new byte[12];
    oct2bin("7777", result, 0, 0);
    assert(result, cast(byte[])[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]);

    result[] = cast(byte)0;

    oct2bin("7776", result[0 .. 11], 0, 0);
    assert(result, cast(byte[])[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1]);
}



void hex2bin(string oct, short[] bin, int skiplast, int flip)
in{
    assert(bin.length <= oct.length * 4);
}
body{
    auto tmp = new short[oct.length * 4];

    {
        auto sink = tmp;
        sink.put(oct.map!"a - '0'"().toBinaryDigits!false(4).map!"cast(short)(a ? -1 : 1)"());
    }

    if(skiplast)
        bin[] = tmp[0 .. bin.length];
    else
        bin[] = tmp[$ - bin.length .. $];

    if(flip)
        bin.reverse();
}
unittest{
    auto result = new byte[12];
    oct2bin("FFFF", result, 0, 0);
    assert(result, cast(byte[])[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]);

    result[] = cast(byte)0;

    oct2bin("FFFE", result[0 .. 11], 0, 0);
    assert(result, cast(byte[])[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1]);
}



/**
渡されたビット列を表すレンジの要素を，nbitまとめたレンジを返します．
packBinaryDigitsとtoBinaryDigitsは互いに逆関数の関係にあります．
*/
auto packBinaryDigits(I, bool fromLSB = true, R)(R rng, size_t nDigits)
if(isInputRange!R && isIntegral!I && is(Unqual!(ElementType!R) == bool))
in{
    assert(nDigits <= I.sizeof * 8);
}
body{
    static struct BinaryDigitsPacker()
    {
        I front() const { return _front; }


        void popFront()
        {
            _front = 0;
            foreach(i; 0 .. _nDigs){
                if(_rng.empty){
                    _empty = true;
                    return;
                }

                if(_rng.front)
                    _front += 1 << i;

                _rng.popFront();
            }
        }


        bool empty() const @property { return _empty; }


      private:
        R _rng;
        I _front;
        size_t _nDigs;
        bool _empty;
    }


    return BinaryDigitsPacker!()(rng, 0, nDigits, false);
}

