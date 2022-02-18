module dffdd.math.exprtemplate;

import std.experimental.allocator;
import std.typecons;
import std.traits;
import std.meta;

import mir.ndslice;

import dffdd.math.vector;
import dffdd.math.matrix;
import dffdd.math.complex;
import dffdd.math.linalg : isBlasType;

/** 
 * 
 */
private
auto
    makeViewOrNewSlice(Vec, Alloc)(const Vec vec, ref Alloc alloc)
if(isVectorLike!Vec)
{
    static if(is(Unqual!(typeof(vec.sliced())) == Slice!(Vec.ElementType*, 1, kind), SliceKind kind))
        return tuple!("view", "isAllocated")(vec.sliced.lightConst, false);
    else
    {
        Slice!(Vec.ElementType*, 1, Contiguous) newslice = alloc.makeArray!(Vec.ElementType)(vec.length).sliced;
        vec.evalTo(newslice, alloc);
        return tuple!("view", "isAllocated")(newslice, true);
    }
}

unittest
{
    auto vec1 = vector!int(3, 1);
    auto s1 = vec1.makeViewOrNewSlice(theAllocator);
    assert(s1.view == [1, 1, 1]);
    assert(s1.isAllocated == false);

    auto vec2 = iota(4).vectored;
    auto s2 = vec2.makeViewOrNewSlice(theAllocator);
    assert(s2.view == [0, 1, 2, 3]);
    assert(s2.isAllocated == true);
    theAllocator.dispose(s2.view.iterator);
}


private
auto
    makeViewOrNewSlice(Mat, Alloc)(const Mat mat, ref Alloc alloc)
if(isMatrixLike!Mat)
{
    static if(is(Unqual!(typeof(mat.sliced())) == Slice!(Mat.ElementType*, 2, kind), SliceKind kind)
      && (kind == Contiguous || kind == Canonical))
    {
        return tuple!("view", "isAllocated")(mat.sliced, false);
    }
    else
    {
        Slice!(Mat.ElementType*, 2, Contiguous) newslice = alloc.makeArray!(Mat.ElementType)(mat.length!0 * mat.length!1).sliced(mat.length!0, mat.length!1);
        mat.evalTo(newslice, alloc);
        return tuple!("view", "isAllocated")(newslice, true);
    }
}



// 
struct MatrixVectorMulGEMV(T, MatA, VecB, VecC)
if(isMatrixLike!MatA && isVectorLike!VecB)
{
    alias ElementType = typeof(T.init * MatA.ElementType.init * VecB.ElementType.init + T.init * VecC.ElementType.init);
    pragma(msg, "E: " ~ ElementType.stringof ~ ", T: " ~ T.stringof ~ ", MatA: " ~ MatA.stringof ~ ", VecB: " ~ VecB.stringof ~ ", VecC: " ~ VecC.stringof);


    this(T alpha, MatA matA, VecB vecB, T beta, VecC vecC)
    {
        _alpha = alpha;
        _beta = beta;
        _matA = matA;
        _vecB = vecB;
        _vecC = vecC;
    }


    ElementType opIndex(size_t n) const
    in(n < this.length)
    {
        ElementType sum = ElementType(0);
        foreach(i; 0 .. _matA.length!1)
            sum += _alpha * _matA[n, i] * _vecB[i];

        return sum + _beta * _vecC[n];
    }


    size_t length() const @property
    {
        return _matA.length!0;
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
    in(dst.length == this.length)
    {
        auto viewA = _matA.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        auto viewB = _vecB.makeViewOrNewSlice(alloc);
        scope(exit) if(viewB.isAllocated) alloc.dispose(viewB.view.iterator);

        auto viewC = _vecC.makeViewOrNewSlice(alloc);
        scope(exit) if(viewC.isAllocated) alloc.dispose(viewC.view.iterator);

        static if(isBlasType!T)
        {
            static if(isFloatingPoint!T || (is(T == MirComplex!E, E) && isFloatingPoint!E))
            {
                dst[] = viewC.view.as!(const(T));

                import mir.blas : gemv;
                gemv(_alpha, viewA.view, viewB.view, _beta, dst);
            }
            else
            {
                dst[] = viewC.view.as!(const(T));
                
                import dffdd.math.linalg : gemv_stdcomplex;
                gemv_stdcomplex(T(_alpha), viewA.view, viewB.view, T(_beta), dst);
            }
        }
        else
        {
            import dffdd.math.linalg : dot;

            foreach(n; 0 .. this.length)
                dst[n] = dot(viewA.view[n].vectored, viewB.view.vectored) * _alpha + _beta * viewC.view[n];
        }
    }


    mixin VectorOperators!(["S*V", "V*S", "V+V"]);


  private:
    T _alpha;
    T _beta;
    MatA _matA;
    VecB _vecB;
    VecC _vecC;
}


MatrixVectorMulGEMV!(T, MatA, VecB, VecC) vectorGemv(T, MatA, VecB, VecC)(T alpha, MatA matA, VecB vecB, T beta, VecC vecC)
{
    return typeof(return)(alpha, matA, vecB, beta, vecC);
}

unittest
{
    static foreach(E; AliasSeq!(/*int, long, float, double, StdComplex!float,*/ MirComplex!double))
    {{
        int err;
        auto mat1 = iota(6).reshape([3, 2], err).as!E.matrixed;
        auto vec1 = [1, 2].sliced.as!E.vectored;
        auto vec2 = [3, 4, 5].sliced.as!E.vectored;
        auto mul1 = vectorGemv(E(2), mat1, vec1, E(3), vec2);

        assert(err == 0);
        assert(mul1.length == 3);
        assert(mul1[0] == 13);
        assert(mul1[1] == 28);
        assert(mul1[2] == 43);

        auto s = slice!E(3);
        mul1.evalTo(s, theAllocator);
        assert(s[0] == 13);
        assert(s[1] == 28);
        assert(s[2] == 43);

        auto view = mul1.makeViewOrNewSlice(theAllocator);
        assert(view.view == [13, 28, 43]);
        assert(view.isAllocated == true);
        static assert(isVectorLike!(typeof(mul1)));
    }}
}


// aX + Y
struct VectorAxpy(T, VecA, VecB)
if(isVectorLike!VecA && isVectorLike!VecB)
{
    alias ElementType = typeof(T.init * VecA.ElementType.init + VecB.ElementType.init);

    this(T alpha, VecA vecA, VecB vecB)
    {
        _alpha = alpha;
        _vecA = vecA;
        _vecB = vecB;
    }


    ElementType opIndex(size_t n) const
    in(n < this.length)
    {
        return _alpha * _vecA[n] + _vecB[n];
    }


    size_t length() const
    {
        return _vecA.length;
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
    in(dst.length == this.length)
    {
        auto viewA = _vecA.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        auto viewB = _vecB.makeViewOrNewSlice(alloc);
        scope(exit) if(viewB.isAllocated) alloc.dispose(viewB.view.iterator);

        dst[] = viewA.view;
        dst[] *= _alpha;
        dst[] += viewB.view;
    }


    mixin VectorOperators!(["S*V", "V*S", "V+V"]);


  private:
    T _alpha;
    VecA _vecA;
    VecB _vecB;
}


VectorAxpy!(T, VecA, VecB) vectorAxpy(T, VecA, VecB)(T scalar, VecA vecA, VecB vecB)
{
    return typeof(return)(scalar, vecA, vecB);
}


unittest
{
    static foreach(E; AliasSeq!(int, long, float, double, StdComplex!float, MirComplex!double))
    {{
        E alpha = E(2);
        auto vec1 = iota(3).as!E.vectored;
        auto vec2 = [3, 2, 3].sliced.as!E.vectored;

        auto axpy1 = vectorAxpy(alpha, vec1, vec2);
        static assert(isVectorLike!(typeof(axpy1)));

        assert(axpy1.length == 3);
        assert(axpy1[0] == 3);
        assert(axpy1[1] == 4);
        assert(axpy1[2] == 7);

        auto s = slice!E(3);
        axpy1.evalTo(s, theAllocator);
        assert(s[0] == 3);
        assert(s[1] == 4);
        assert(s[2] == 7);
    }}
}



mixin template VectorOperators(string[] list)
{
    import std.algorithm : canFind;
    import std.traits : isFloatingPoint, isIntegral;
    import dffdd.math.linalg : isBlasType;

    static if(canFind(list, "S*V"))
    {
        auto opBinary(string op : "*", T)(T scalar)
        if(isIntegral!T || isFloatingPoint!T || isBlasType!T)
        {
            return vectorAxpy(scalar, this, 0.repeat(this.length).vectored);
        }
    }


    static if(canFind(list, "V*S"))
    {
        auto opBinaryRight(string op : "*", T)(T scalar)
        if(isIntegral!T || isFloatingPoint!T || isBlasType!T)
        {
            return vectorAxpy(scalar, this, 0.repeat(this.length).vectored);
        }
    }


    static if(canFind(list, "V+V"))
    {
        auto opBinary(string op : "+", V)(V vec)
        if(isVectorLike!V)
        in(vec.length == this.length)
        {
            return vectorAxpy(ElementType(1), this, vec);
        }
    }
}


mixin template MatrixOperators(string[] list)
{
    static if(canFind(list, "M*V"))
    {
        auto opBinary(string op : "*", V)(V vec)
        if(isVectorLike!V)
        in(vec.length == this.length!1)
        {
            return vectorGemv(ElementType(1), this, vec, ElementType(0), 0.repeat(this.length!0).vectored);
        }
    }
}
