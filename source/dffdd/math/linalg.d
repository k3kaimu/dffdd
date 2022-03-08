module dffdd.math.linalg;

import std.traits;
import std.experimental.allocator;

import dffdd.math.complex;
import dffdd.math.vector;
import dffdd.math.matrix;
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
        alias Ret = Unqual!(typeof(vecA[0] * vecB[0]));

        static if(isComplex!Ret)
            return reduce!((a, b, c) => a + conj(b) * c)(Ret(0), vecA.sliced.as!Ret, vecB.sliced.as!Ret);
        else
            return reduce!((a, b, c) => a + b * c)(Ret(0), vecA.sliced.as!Ret, vecB.sliced.as!Ret);
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


/**
複素行列を以下のような実数行列へ変換します
[[R  -I],
 [I   R]]
ここで，RはRe{M}でIはIm{M}です
*/
RealMatrixRIIRStructure!Mat toRealRIIR(Mat)(Mat mat)
if(isMatrixLike!Mat && isComplex!(Mat.ElementType))
{
    return RealMatrixRIIRStructure!Mat(mat);
}


struct RealMatrixRIIRStructure(Mat)
if(isMatrixLike!Mat && isComplex!(Mat.ElementType))
{
    alias ElementType = typeof(Mat.ElementType.init.re);
    enum exprTreeCost = Mat.exprTreeCost + 1;


    this(Mat mat)
    {
        _mat = mat;
    }


    size_t length(size_t dim)() const {
        return _mat.length!dim * 2;
    }


    ElementType opIndex(size_t i, size_t j) const
    {
        if(i < _mat.length!0 && j < _mat.length!1)
            return _mat[i, j].re;
        else if(i < _mat.length!0 && j >= _mat.length!1)
            return -_mat[i, j - _mat.length!1].im;
        else if(i >= _mat.length!0 && j < _mat.length!1)
            return _mat[i - _mat.length!0, j].im;
        else
            return _mat[i - _mat.length!0, j - _mat.length!1].re;
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        auto viewA = _mat.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        auto v = viewA.view.lightConst;
        immutable N = v.length!0,
                  M = v.length!1;

        dst[0 .. N  , 0 .. M  ] = v.map!"a.re";
        dst[0 .. N  , M .. M*2] = v.map!"-a.im";
        dst[N .. N*2, 0 .. M  ] = v.map!"a.im";
        dst[N .. N*2, M .. M*2] = v.map!"a.re";
    }


    import dffdd.math.exprtemplate;
    mixin MatrixOperators!(["M*M", "M*V"]) OpImpls;
    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*S", ".T", ".H"]));


    auto _opB_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(this.length!1 == mat.length!0)
    {
        static if(is(typeof(mat.fromRealRIIR())))
        {
            alias T1 = typeof(toRealRIIR(_mat * mat.fromRealRIIR));
            alias T2 = typeof(OpImpls._opB_Impl_!"*"(mat));

            static if(T1.exprTreeCost < T2.exprTreeCost)
                return toRealRIIR(_mat * mat.fromRealRIIR);
            else
                return OpImpls._opB_Impl_!"*"(mat);
        }
        else
            return OpImpls._opB_Impl_!"*"(mat);
    }


    auto _opB_Impl_(string op : "*", V)(V vec)
    if(isVectorLike!V)
    in(this.length!1 == vec.length)
    {
        static if(is(typeof(vec.fromRealRI())))
        {
            alias T1 = typeof(toRealRI(_mat * vec.fromRealRI));
            alias T2 = typeof(OpImpls._opB_Impl_!"*"(vec));

            static if(T1.exprTreeCost < T2.exprTreeCost)
                return toRealRI(_mat * vec.fromRealRI);
            else
                return OpImpls._opB_Impl_!"*"(vec);
        }
        else
            return OpImpls._opB_Impl_!"*"(vec);
    }


    Mat fromRealRIIR()
    {
        return _mat;
    }
    

  private:
    Mat _mat;
}

unittest
{
    alias C = MirComplex!float;
    auto cmat = [C(1, 2)].sliced(1, 1).matrixed;
    auto rmat = cmat.toRealRIIR;
    assert(rmat[0, 0] == 1);
    assert(rmat[0, 1] == -2);
    assert(rmat[1, 0] == 2);
    assert(rmat[1, 1] == 1);

    auto s = slice!float(2, 2);
    rmat.evalTo(s, matvecAllocator);
    assert(s == [[1, -2], [2, 1]]);
}

unittest
{
    alias C = MirComplex!float;
    auto cmat = [C(1, 2)].sliced(1, 1).matrixed;
    auto rmat = cmat.toRealRIIR;
    
    auto cvec = [C(1, 2)].sliced(1).vectored;
    auto rvec = cvec.toRealRI;

    auto rmul = rmat * rvec;
    auto cmul = rmul.fromRealRI;

    import dffdd.math.exprtemplate;
    static assert(is(typeof(cmul) == MatrixVectorMulGEMV!(C, M, V1, V2), M, V1, V2));


    auto rmul2 = rmat * rmat;
    auto cmul2 = rmul2.fromRealRIIR;
    static assert(is(typeof(cmul2) == MatrixMatrixMulGEMM!(C, M1, M2, M3), M1, M2, M3));
}



/**
複素ベクトルを[R I]という実数ベクトルへ変換します．
*/
auto toRealRI(Vec)(Vec vec)
if(isVectorLike!Vec && isComplex!(Vec.ElementType))
{
    static if(is(Vec == ComplexVectorFromRIStructure!M, M))
        return vec._vec;
    else
        return RealVectorRIStructure!Vec(vec);
}


struct RealVectorRIStructure(Vec)
if(isVectorLike!Vec && isComplex!(Vec.ElementType))
{
    alias ElementType = typeof(Vec.ElementType.init.re);
    enum exprTreeCost = Vec.exprTreeCost + 1;


    this(Vec vec)
    {
        _vec = vec;
    }


    size_t length() const {
        return _vec.length * 2;
    }


    ElementType opIndex(size_t i) const
    {
        if(i < _vec.length)
            return _vec[i].re;
        else
            return _vec[i - _vec.length].im;
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
    in(dst.length == this.length)
    {
        auto viewA = _vec.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        auto v = viewA.view.lightConst;
        immutable N = v.length;

        dst[0 .. N  ] = v.map!"a.re";
        dst[N .. N*2] = v.map!"a.im";
    }


    import dffdd.math.exprtemplate;
    mixin(definitionsOfVectorOperators(["defaults", "V+V", "V*S"]));

  private:
    Vec _vec;
}


unittest
{
    alias C = MirComplex!float;
    auto cvec = [C(1, 2)].sliced(1).vectored;
    auto rvec = cvec.toRealRI;
    assert(rvec[0] == 1);
    assert(rvec[1] == 2);

    auto s = slice!float(2);
    rvec.evalTo(s, matvecAllocator);
    assert(s == [1, 2]);
}



/**
実数ベクトルを複素ベクトルへ変換します
*/
auto fromRealRI(alias C = MirComplex, Vec)(Vec vec)
if(isVectorLike!Vec && !isComplex!(Vec.ElementType))
in(vec.length % 2 == 0)
{
    static if(is(Vec == RealVectorRIStructure!V, V) && is(V.ElementType == C!(Vec.ElementType)))
        return vec._vec;
    else
        return ComplexVectorFromRIStructure!(C!(Vec.ElementType), Vec)(vec);
}


struct ComplexVectorFromRIStructure(C, Vec)
if(isVectorLike!Vec && !isComplex!(Vec.ElementType) && isComplex!C)
{
    alias ElementType = C;
    enum exprTreeCost = Vec.exprTreeCost + 1;


    this(Vec vec)
    {
        _vec = vec;
    }


    size_t length() const {
        return _vec.length / 2;
    }


    ElementType opIndex(size_t i) const
    {
        return C(_vec[i], _vec[i + _vec.length/2]);
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
    in(dst.length == this.length)
    {
        auto viewA = _vec.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        auto v = viewA.view.lightConst;
        immutable N = v.length;

        dst[] += v[0 .. $/2].map!(a => T(a, 0));
        dst[] += v[$/2 .. $].map!(a => T(0, a));
    }


    import dffdd.math.exprtemplate;
    mixin(definitionsOfVectorOperators(["defaults", "V+V", "V*S"]));

  private:
    Vec _vec;
}

unittest
{
    alias C = MirComplex!float;
    auto rvec = [1, 2].sliced(2).as!float.vectored;
    assert(rvec[0] == 1);
    assert(rvec[1] == 2);

    auto cvec = rvec.fromRealRI;
    assert(cvec[0] == C(1, 2));
    auto s = slice!C(1);
    cvec.evalTo(s, matvecAllocator);
    assert(s == [C(1, 2)]);
}

unittest
{
    alias C = MirComplex!float;
    auto cvec1 = [C(1, 2)].sliced(1).vectored;
    auto cvec2 = cvec1.toRealRI.fromRealRI;
    static assert(is(typeof(cvec2) == typeof(cvec1)));
}


unittest
{
    alias C = MirComplex!float;
    auto cmat = [C(1, 2)].sliced(1, 1).matrixed;
    auto rmat = cmat.toRealRIIR;
    assert(rmat[0, 0] == 1);
    assert(rmat[0, 1] == -2);
    assert(rmat[1, 0] == 2);
    assert(rmat[1, 1] == 1);

    auto cvec = [C(1, 2)].sliced(1).vectored;
    auto rvec = cvec.toRealRI;
    assert(rvec[0] == 1);
    assert(rvec[1] == 2);

    auto _1 = rmat._opB_Impl_!"*"(rvec);
    auto rmul = rmat * rvec;
    auto cmul = rmul.fromRealRI;
    assert(cmul[0] == C(1, 2) * C(1, 2));

    auto s = slice!C(1);
    cmul.evalTo(s, matvecAllocator);
    assert(s == [C(1, 2) * C(1, 2)]);
}


void qrDecomp(MatA, C, SliceKind kindQ, SliceKind kindR, Alloc)(in MatA matA, MatrixedSlice!(C*, kindQ) matQ, MatrixedSlice!(C*, kindR) matR, ref Alloc alloc = matvecAllocator)
in(matA.length!0 == matQ.length!0)
in(matQ.length!0 == matQ.length!1)
in(matA.length!0 == matR.length!0)
in(matA.length!1 == matR.length!1)
{
    import kaleidic.lubeck : qrDecomp;

    auto viewA = makeViewOrNewSlice(matA, alloc);
    scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

    auto qrResult = qrDecomp(viewA.view.as!C.slice);
    matQ.sliced()[] = qrResult.Q;
    matR.sliced()[] = qrResult.R;
}

unittest
{
    import std.math : isClose;

    alias C = MirComplex!double;
    auto mat = matrix!C(3, 3);
    mat.sliced()[] = [[C(12), C(-51),   C(4)],
              [C(6), C(167), C(-68)],
              [C(-4),  C(24), C(-41)]];

    auto matQ = matrix!C(3, 3);
    auto matR = matrix!C(3, 3);

    mat.qrDecomp(matQ, matR);
    auto mul = matrix!C(3, 3);
    mul[] = matQ * matR;
    foreach(i; 0 .. 3)
        foreach(j; 0 .. 3) {
            assert(mul[i, j].re.isClose(mat[i, j].re));
            assert(mul[i, j].im.isClose(mat[i, j].im));
        }
}


struct InverseMatrix(M)
if(isMatrixLike!M)
{
    alias ElementType = M.ElementType;
    enum exprTreeCost = M.exprTreeCost + EXPR_COST_N^^3;


    this(M mat)
    {
        _mat = mat;
    }


    size_t length(size_t dim)() const { return _mat.length!dim; }


    ElementType opIndex(size_t, size_t) const { assert(0); return ElementType(float.nan); }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        _mat.evalTo(dst, alloc);
        invImpl(dst, alloc);
    }


    // import dffdd.math.exprtemplate;
    // mixin(definitionsOfVectorOperators(["defaults", "V+V", "V*S"]));


  private:
    M _mat;
}


auto inv(M)(M mat)
{
    static if(is(typeof(mat._inv_Impl_())))
        return mat._inv_Impl_;
    else
        return InverseMatrix!M(mat);
}


/**
逆行列
*/
void invImpl(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) matA, ref Alloc alloc)
if(isBlasType!T)
in{
    assert(matA.length!0 == matA.length!1);
}
do{
    static if(kindA == Universal)
    {
        if (matA._stride!1 != 1)
        {
            assert(matA._stride!0 == 1, "Matrix A must have a stride equal to 1.");
            .invImpl(
                matA.universal.transposed,
                alloc);
            return;
        }
        assert(matA._stride!1 == 1, "Matrix A must have a stride equal to 1.");
    }

    int[] ipiv = alloc.makeArray!int(matA.length!0);
    scope(exit) alloc.dispose(ipiv);

    immutable n = matA.length!0;
    immutable lda = matA.matrixStride;

    import dffdd.utils.linalg : LAPACKE_sgetrf, LAPACKE_sgetri, LAPACKE_dgetrf, LAPACKE_dgetri,
                                LAPACKE_cgetrf, LAPACKE_cgetri, LAPACKE_zgetrf, LAPACKE_zgetri;

    import cblas : Order;

    static if(is(T == float))
    {
        alias getrf = LAPACKE_sgetrf;
        alias getri = LAPACKE_sgetri;
        alias E = float;
    }
    else static if(is(T == double))
    {
        alias getrf = LAPACKE_dgetrf;
        alias getri = LAPACKE_dgetri;
        alias E = double;
    }
    else static if(is(typeof(T.init.re) == float))
    {
        alias getrf = LAPACKE_cgetrf;
        alias getri = LAPACKE_cgetri;
        alias E = float[2];
    }
    else static if(is(typeof(T.init.re) == double))
    {
        alias getrf = LAPACKE_zgetrf;
        alias getri = LAPACKE_zgetri;
        alias E = double[2];
    }
    else static assert(0);

    getrf(Order.RowMajor, cast(uint)n, cast(uint)n, cast(E*)&(matA[0, 0]), cast(uint)lda, ipiv.ptr);
    getri(Order.RowMajor, cast(uint)n, cast(E*)&(matA[0, 0]), cast(uint)lda, ipiv.ptr);
}

unittest
{
    import std.math;
    import std.stdio;
    auto mat = new Complex!float[4].sliced(2, 2);
    mat[0, 0] = 1;
    mat[0, 1] = 0;
    mat[1, 0] = 0;
    mat[1, 1] = 1;

    invImpl(mat, matvecAllocator);
    assert(isClose(mat[0, 0].re, 1));
    assert(isClose(mat[0, 1].re, 0));
    assert(isClose(mat[1, 0].re, 0));
    assert(isClose(mat[1, 1].re, 1));
    foreach(i; [0, 1]) foreach(j; [0, 1]) assert(isClose(mat[i, j].im, 0));


    mat[0, 0] = 3;
    mat[0, 1] = 4;
    mat[1, 0] = 1;
    mat[1, 1] = 2;
    invImpl(mat, matvecAllocator);
    //writeln(mat);

    assert(isClose(mat[0, 0].re, 1));
    assert(isClose(mat[0, 1].re, -2));
    assert(isClose(mat[1, 0].re, -0.5));
    assert(isClose(mat[1, 1].re, 1.5));
    foreach(i; [0, 1]) foreach(j; [0, 1]) assert(isClose(mat[i, j].im, 0));
}


unittest
{
    import std.math;
    import std.stdio;
    auto mat = new float[4].sliced(2, 2);
    mat[0, 0] = 1;
    mat[0, 1] = 0;
    mat[1, 0] = 0;
    mat[1, 1] = 1;

    invImpl(mat, matvecAllocator);
    assert(isClose(mat[0, 0], 1));
    assert(isClose(mat[0, 1], 0));
    assert(isClose(mat[1, 0], 0));
    assert(isClose(mat[1, 1], 1));


    mat[0, 0] = 3;
    mat[0, 1] = 4;
    mat[1, 0] = 1;
    mat[1, 1] = 2;
    invImpl(mat, matvecAllocator);
    //writeln(mat);

    assert(isClose(mat[0, 0], 1));
    assert(isClose(mat[0, 1], -2));
    assert(isClose(mat[1, 0], -0.5));
    assert(isClose(mat[1, 1], 1.5));
}


unittest
{
    import std.math;
    import std.stdio;
    auto mat = new float[4].sliced(2, 2).matrixed;
    mat[0, 0] = 1;
    mat[0, 1] = 0;
    mat[1, 0] = 0;
    mat[1, 1] = 1;

    // invImpl(mat, matvecAllocator);
    // auto inv = InverseMatrix!(typeof(mat))(mat);
    auto inv = mat.inv;
    assert(isClose(mat[0, 0], 1));
    assert(isClose(mat[0, 1], 0));
    assert(isClose(mat[1, 0], 0));
    assert(isClose(mat[1, 1], 1));


    mat[0, 0] = 3;
    mat[0, 1] = 4;
    mat[1, 0] = 1;
    mat[1, 1] = 2;
    invImpl(mat, matvecAllocator);
    //writeln(mat);

    assert(isClose(mat[0, 0], 1));
    assert(isClose(mat[0, 1], -2));
    assert(isClose(mat[1, 0], -0.5));
    assert(isClose(mat[1, 1], 1.5));
}


