module dffdd.math.linalg;

import std.traits;
import std.experimental.allocator;

import dffdd.math.complex;
import dffdd.math.vector;
import dffdd.math.exprtemplate : makeViewOrNewSlice, matvecAllocator;

import mir.ndslice;


// copy from https://github.com/libmir/mir-blas/blob/master/source/mir/blas.d
private auto matrixStride(S)(S a)
 if (S.N == 2)
{
    assert(a._stride!1 == 1);
    return a._stride != 1 ? a._stride : a.length!1;
}


enum isBlasType(T)
                 = is(Unqual!T == float) || is(Unqual!T == double)
                || is(Unqual!T == StdComplex!float) || is(Unqual!T == StdComplex!double)
                || is(Unqual!T == MirComplex!float) || is(Unqual!T == MirComplex!double);

unittest
{
    static assert(isBlasType!float);
    static assert(isBlasType!(StdComplex!float));
    static assert(isBlasType!(MirComplex!float));
}


auto dotImpl(IteratorA, SliceKind kindA, IteratorB, SliceKind kindB)
    (VectoredSlice!(IteratorA, kindA) vecA, VectoredSlice!(IteratorB, kindB) vecB)
in(vecA.length == vecB.length)
{
    static if(is(IteratorA == E*, E) && is(IteratorB == E*) && isBlasType!E)
    {
        import cblas;

        auto sliceA = vecA.sliced;
        auto sliceB = vecB.sliced;

        static if(is(Unqual!E == float))
            return cblas_sdot(cast(blasint)sliceA.length, sliceA.iterator, cast(blasint)sliceA._stride!0, sliceB.iterator, cast(blasint)sliceB._stride!0);
        else static if(is(Unqual!E == double))
            return cblas_ddot(cast(blasint)sliceA.length, sliceA.iterator, cast(blasint)sliceA._stride!0, sliceB.iterator, cast(blasint)sliceB._stride!0);
        else static if(is(Unqual!(typeof(E.init.re)) == float))
        {
            auto ret = cblas_cdotu(cast(blasint)sliceA.length,
                cast(MirComplex!float*)sliceA.iterator, cast(blasint)sliceA._stride!0,
                cast(MirComplex!float*)sliceB.iterator, cast(blasint)sliceB._stride!0);

            static if(is(E == cfloat))
                return ret.re + ret.im*1.0fi;
            else
                return E(ret.re, ret.im);
        }
        else static if(is(Unqual!(typeof(E.init.re)) == double))
        {
            auto ret = cblas_zdotu(cast(blasint)sliceA.length,
                cast(MirComplex!double*)sliceA.iterator, cast(blasint)sliceA._stride!0,
                cast(MirComplex!double*)sliceB.iterator, cast(blasint)sliceB._stride!0);

            return E(ret.re, ret.im);
        }
        else
            static assert(0, E.stringof ~ " is invalid type for dot product.");
    }
    else
    {
        alias Ret = typeof(vecA[0] * vecB[0]);
        return reduce!"a + b * c"(Ret(0), vecA.sliced.as!Ret, vecB.sliced.as!Ret);
    }
}

unittest
{
    auto vec1 = iota([3], 0, 1).map!"a*2".as!float.universal.vector;
    assert(dotImpl(vec1, vec1) == 2*2+4*4);

    auto vec2 = iota([3], 0, 1).map!(a => MirComplex!float(a*2, 0)).vector;
    assert(dotImpl(vec2, vec2) == 2*2+4*4);

    auto vec3 = [MirComplex!float(1, 1), MirComplex!float(1, -2)].sliced.vectored;
    assert(dotImpl(vec3, vec3) == MirComplex!float(1, 1) * MirComplex!float(1, 1) + MirComplex!float(1, -2) * MirComplex!float(1, -2));
}


auto dot(VecA, VecB, Alloc)(VecA vecA, VecB vecB, ref Alloc alloc = matvecAllocator)
if(isVectorLike!VecA && isVectorLike!VecB)
in(vecA.length == vecB.length)
{
    auto viewA = vecA.makeViewOrNewSlice(alloc);
    scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

    auto viewB = vecB.makeViewOrNewSlice(alloc);
    scope(exit) if(viewB.isAllocated) alloc.dispose(viewB.view.iterator);

    return dotImpl(viewA.view.vectored, viewB.view.vectored);
}



auto dotHImpl(IteratorA, SliceKind kindA, IteratorB, SliceKind kindB)
    (VectoredSlice!(IteratorA, kindA) vecA, VectoredSlice!(IteratorB, kindB) vecB)
in(vecA.length == vecB.length)
{
    static if(is(IteratorA == E*, E) && is(IteratorB == E*) && isBlasType!E)
    {
        import cblas;

        static if(is(Unqual!E == float) || is(Unqual!E == double))
            return .dot(vecA, vecB);
        else {
            auto sliceA = vecA.sliced;
            auto sliceB = vecB.sliced;

            static if(is(Unqual!(typeof(E.init.re)) == float))
            {
                auto ret = cblas_cdotc(cast(blasint)sliceA.length,
                    cast(MirComplex!float*)sliceA.iterator, cast(blasint)sliceA._stride!0,
                    cast(MirComplex!float*)sliceB.iterator, cast(blasint)sliceB._stride!0);

                static if(is(Unqual!E == cfloat))
                    return ret.re + ret.im*1.0fi;
                else
                    return E(ret.re, ret.im);
            }
            else static if(is(Unqual!(typeof(E.init.re)) == double))
            {
                auto ret = cblas_zdotc(cast(blasint)sliceA.length,
                    cast(MirComplex!double*)sliceA.iterator, cast(blasint)sliceA._stride!0,
                    cast(MirComplex!double*)sliceB.iterator, cast(blasint)sliceB._stride!0);

                return E(ret.re, ret.im);
            }
            else
                static assert(0, E.stringof ~ " is invalid type for dot product.");
        }
    }
    else
    {
        alias Ret = typeof(vecA[0] * vecB[0]);
        return reduce!(a => a + conj(b) * c)(Ret(0), vecA.sliced.as!Ret, vecB.sliced.as!Ret);
    }
}


auto dotH(VecA, VecB, Alloc)(VecA vecA, VecB vecB, ref Alloc alloc = matvecAllocator)
if(isVectorLike!VecA && isVectorLike!VecB)
in(vecA.length == vecB.length)
{
    auto viewA = vecA.makeViewOrNewSlice(alloc);
    scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

    auto viewB = vecB.makeViewOrNewSlice(alloc);
    scope(exit) if(viewB.isAllocated) alloc.dispose(viewB.view.iterator);

    return dotHImpl(viewA.view.vectored, viewB.view.vectored);
}


unittest
{
    auto vec1 = iota([3], 0, 1).map!"a*2".as!float.universal.vector;
    assert(dotHImpl(vec1, vec1) == 2*2+4*4);

    auto vec2 = iota([3], 0, 1).map!(a => MirComplex!float(a*2, 0)).vector;
    assert(dotHImpl(vec2, vec2) == 2*2+4*4);

    auto vec3 = [MirComplex!float(1, 1), MirComplex!float(1, -2)].sliced.vectored;
    assert(dotHImpl(vec3, vec3) == MirComplex!float(1, 1) * MirComplex!float(1, -1) + MirComplex!float(1, -2) * MirComplex!float(1, 2));
}


auto norm2Impl(IteratorA, SliceKind kindA)(VectoredSlice!(IteratorA, kindA) vec)
{
    return dotHImpl(vec, vec);
}


auto norm2(VecA, Alloc)(VecA vecA, ref Alloc alloc = matvecAllocator)
if(isVectorLike!VecA)
in(vecA.length)
{
    auto viewA = vecA.makeViewOrNewSlice(alloc);
    scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

    return norm2Impl(viewA.view.vectored);
}


unittest
{
    auto vec3 = [MirComplex!float(1, 1), MirComplex!float(1, -2)].sliced.vectored;
    assert(norm2(vec3) == 1+1+1+4);
}


auto meanImpl(Iterator, SliceKind kind)(VectoredSlice!(Iterator, kind) vec)
if(isVectorLike!VecA)
{
    alias E = vec.ElementType;
    return reduce!"a + b"(E(0), vec.sliced) / vec.length;
}


auto mean(VecA)(VecA vecA, ref Alloc alloc = matvecAllocator)
if(isVectorLike!VecA)
{
    auto viewA = vecA.makeViewOrNewSlice(alloc);
    scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

    return meanImpl(viewA.view.vectored);
}


// copy and modified from https://github.com/libmir/mir-blas/blob/master/source/mir/blas.d
void gemv_stdcomplex(T, SliceKind kindA, SliceKind kindX, SliceKind kindY)
    (T alpha, Slice!(const(T)*, 2, kindA) a, Slice!(const(T)*, 1, kindX) x, T beta, Slice!(T*, 1, kindY) y)
if(is(T == StdComplex!E, E))
in
{
    assert(a.length!1 == x.length);
    assert(a.length!0 == y.length);
}
do
{
    alias E = typeof(T.init.re);

    static if (kindA == Universal)
    {
        bool transA;
        if (a._stride!1 != 1)
        {
            a = a.transposed;
            transA = true;
        }
        assert(a._stride!1 == 1, "Matrix A must have a stride equal to 1.");
    }
    else
        enum transA = false;

    import cblas;
    immutable
        _alpha = MirComplex!E(alpha.re, alpha.im),
        _beta = MirComplex!E(beta.re, beta.im);

    gemv(
        transA ? cblas.Order.ColMajor : cblas.Order.RowMajor,
        cblas.Transpose.NoTrans,
        
        cast(cblas.blasint) y.length,
        cast(cblas.blasint) x.length,

        _alpha,

        cast(const(MirComplex!E)*)a.iterator,
        cast(cblas.blasint) a.matrixStride,

        cast(const(MirComplex!E)*)x.iterator,
        cast(cblas.blasint) x._stride,

        _beta,

        cast(MirComplex!E*)y.iterator,
        cast(cblas.blasint) y._stride,
    );
}


// copy and modified from https://github.com/libmir/mir-blas/blob/master/source/mir/blas.d
void gemm_stdcomplex(T, SliceKind kindA, SliceKind kindB, SliceKind kindC,)
    (T alpha, Slice!(const(T)*, 2, kindA) a, Slice!(const(T)*, 2, kindB) b, T beta, Slice!(T*, 2, kindC) c)
if(is(T == StdComplex!E, E))
in
{
    assert(a.length!1 == b.length!0);
    assert(a.length!0 == c.length!0);
    assert(c.length!1 == b.length!1);
}
do
{
    import cblas;
    alias E = typeof(T.init.re);

    auto k = cast(cblas.blasint) a.length!1;

    static if (kindC == Universal)
    {
        if (c._stride!1 != 1)
        {
            assert(c._stride!0 == 1, "Matrix C must have a stride equal to 1.");
            .gemm(
                alpha,
                b.universal.transposed,
                a.universal.transposed,
                beta,
                c.transposed.assumeCanonical);
            return;
        }
        assert(c._stride!1 == 1, "Matrix C must have a stride equal to 1.");
    }
    static if (kindA == Universal)
    {
        bool transA;
        if (a._stride!1 != 1)
        {
            a = a.transposed;
            transA = true;
        }
        assert(a._stride!1 == 1, "Matrix A must have a stride equal to 1.");
    }
    else
        enum transA = false;
    static if (kindB == Universal)
    {
        bool transB;
        if (b._stride!1 != 1)
        {
            b = b.transposed;
            transB = true;
        }
        assert(b._stride!1 == 1, "Matrix B must have a stride equal to 1.");
    }
    else
        enum transB = false;

    import cblas;
    immutable
        _alpha = MirComplex!E(alpha.re, alpha.im),
        _beta = MirComplex!E(beta.re, beta.im);

    gemm(
        cblas.Order.RowMajor,
        transA ? cblas.Transpose.Trans : cblas.Transpose.NoTrans,
        transB ? cblas.Transpose.Trans : cblas.Transpose.NoTrans,
        
        cast(cblas.blasint) c.length!0,
        cast(cblas.blasint) c.length!1, 
        cast(cblas.blasint) k,

        _alpha,

        cast(const(MirComplex!E)*)a.iterator,
        cast(cblas.blasint) a.matrixStride,

        cast(const(MirComplex!E)*)b.iterator,
        cast(cblas.blasint) b.matrixStride,

        _beta,

        cast(MirComplex!E*)c.iterator,
        cast(cblas.blasint) c.matrixStride,
    );
}