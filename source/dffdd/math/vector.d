module dffdd.math.vector;

import std.traits;

import mir.ndslice;
import dffdd.math.linalg;
import std.experimental.allocator;
import dffdd.math.exprtemplate : isExpressionTemplate;



VectoredSlice!(Iterator, kind) vectored(Iterator, SliceKind kind)(Slice!(Iterator, 1, kind) slice)
{
    return VectoredSlice!(Iterator, kind)(slice);
}


Vector!(T, Contiguous) vector(T)(size_t len)
{
    return vector!T(len, T.init);
}


Vector!(T, Contiguous) vector(T)(size_t len, T init)
{
    return Vector!(T, Contiguous)(slice!T([len], init));
}


auto vector(Iterator, SliceKind kind)(Slice!(Iterator, 1, kind) s)
{
    return vectored(slice(s));
}


auto makeVector(T, Alloc)(ref Alloc alloc, size_t size)
{
    T[] arr = alloc.makeArray!(T)(size);
    return arr.sliced.vectored;
}


auto makeVector(T, Alloc)(ref Alloc alloc, size_t size, T init)
{
    T[] arr = alloc.makeArray!(T)(size);
    arr[] = init;
    return arr.sliced.vectored;
}


enum isVectorLike(T) = isExpressionTemplate!T && is(typeof((T t, size_t i){
    T.ElementType e = t[i];
    size_t len = t.length;

    auto dst = slice!(typeof(e))(len);
    auto allocator = std.experimental.allocator.theAllocator;
    t.evalTo(dst, allocator);
}));


struct VectoredSlice(Iterator, SliceKind kind)
{
    alias ElementType = Unqual!(typeof(Slice!(Iterator, 1, kind).init[0]));
    enum float exprTreeDepth = 0;


    this(Slice!(Iterator, 1, kind) s)
    {
        _slice = s;

        static if(is(Iterator == ElementType*))
        {
            _tmp = slice!ElementType(s.length);
        }
    }


    size_t length() const @property { return _slice.length; }


    Slice!(Iterator, 1, kind) sliced() { return _slice; }
    auto sliced() const { return _slice.lightConst; }


    auto lightConst() const
    {
        return this.sliced.lightConst.vectored;
    }


    auto ref opIndex(size_t n)
    {
        return _slice[n];
    }


    auto opIndex(size_t n) const
    {
        return _slice[n];
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
    in(dst.length == this.length)
    {
        dst[] = _slice.lightConst.as!(const(T));
    }


    import dffdd.math.exprtemplate;
    mixin(definitionsOfVectorOperators(["defaults", "V*S", "V+V"]));

    static if(is(Iterator == ElementType*))
    {
        auto opSliceAssign(V)(V vec)
        if(isVectorLike!V)
        in(vec.length == this.length)
        {
            vec.evalTo(_tmp, matvecAllocator);
            _slice[] = _tmp;
        }


        auto opSliceOpAssign(string op, V)(V vec)
        if(isVectorLike!V && (op == "+" || op == "-"))
        in(vec.length == this.length)
        {
            vec.evalTo(_tmp, matvecAllocator);

            static if(op == "+")
                _slice[] += _tmp;
            else
                _slice[] -= _tmp;
        }


        auto opSliceAssign(S)(S scalar)
        if(!isMatrixLike!S && !isVectorLike!S)
        {
            this._slice[] = scalar;
        }


        auto opSliceOpAssign(string op, S)(S scalar)
        if(!isMatrixLike!S && !isVectorLike!S && (op == "+" || op == "-"))
        {
            static if(op == "+")
                this._slice[] += scalar;
            else
                this._slice[] -= scalar;
        }
    }


  private:
    Slice!(Iterator, 1, kind) _slice;

    static if(is(Iterator == ElementType*))
    {
        Slice!(ElementType*, 1, Contiguous) _tmp;
    }
}


alias Vector(T, SliceKind kind) = VectoredSlice!(T*, kind);


unittest
{
    auto vec1 = vector!int(4, 1);
    assert(vec1.length == 4);
    assert(vec1[0] == 1);
    assert(vec1[1] == 1);
    assert(vec1[2] == 1);
    assert(vec1[3] == 1);

    vec1[0] = -1;
    assert(vec1[0] == -1);

    auto vec2 = iota([3], 0, 1).map!"a*2".vector;
    assert(vec2[0] == 0);
    assert(vec2[1] == 2);
    assert(vec2[2] == 4);

    static assert(isVectorLike!(typeof(vec1)));
}


unittest
{
    auto vec1 = vector!int(3, 1);
    auto mul1 = vec1 * 2;
    assert(mul1[0] == 2);
    assert(mul1[1] == 2);
    assert(mul1[2] == 2);

    auto vec2 = iota([3], 0, 1).map!"a*2".as!int.vector;
    auto add2 = mul1 + vec2;
    assert(add2[0] == 2);
    assert(add2[1] == 4);
    assert(add2[2] == 6);

    vec1[] = add2;
    assert(vec1[0] == 2);
    assert(vec1[1] == 4);
    assert(vec1[2] == 6);
}
