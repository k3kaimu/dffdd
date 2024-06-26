module dffdd.math.exprtemplate;

import std.algorithm : max, min, sum;
import std.experimental.allocator;
import std.typecons;
import std.traits;
import std.meta;

import mir.ndslice;

import dffdd.math.vector;
import dffdd.math.matrix;
import dffdd.math.matrixspecial;
import dffdd.math.complex;
import dffdd.math.linalg : isBlasType;


alias matvecAllocator = theAllocator;


enum float EXPR_COST_N = 10;


enum isExpressionTemplate(T) = is(typeof(T.exprTreeCost) : float);
enum exprTreeCost(T) = T.exprTreeCost;



class TempMemoryManagerImpl(size_t alignBytes)
{
    import core.memory;
    enum size_t alignment = alignBytes;


    this(size_t initSize)
    {
        _memoryTape = _alignedRealloc(null, initSize, alignment);
        _used = 0;
    }


    CheckPoint saveState() pure nothrow @safe
    {
        return CheckPoint(this, _used);
    }


    E[] makeArray(E)(size_t n) pure nothrow @trusted
    if(!hasIndirections!E)
    {
        size_t last = _used + E.sizeof * n;
        if(_memoryTape.length < last) {
            _memoryTape = _alignedRealloc(_memoryTape, last*2, alignment);
        }

        auto ret = _memoryTape[_used ..last];
        _used = _toNextAlignment(last);
        return cast(E[])ret;
    }


    Slice!(E*, N, Contiguous) makeSlice(E, size_t N)(size_t[N] shape...) pure nothrow @trusted
    if(!hasIndirections!E)
    {
        size_t len = reduce!"a*b"(cast(size_t)1, shape[]);
        auto buf = this.makeArray!E(len);
        return buf.sliced(shape);
    }


    void dispose(E)(E*) pure nothrow @safe {}
    void dispose(E)(E[]) pure nothrow @safe {}


  private:
    void[] _memoryTape;
    size_t _used;


    static struct CheckPoint
    {
        TempMemoryManagerImpl target;
        size_t used;

        ~this()
        {
            target._used = used;
        }
    }


    static size_t _toNextAlignment(size_t n) pure nothrow @safe
    {
        return (n + (alignment - 1)) & ~(alignment - 1);
    }


    static void[] _alignedRealloc(void[] buf, size_t n, size_t alignment) pure nothrow @trusted
    {
        auto newbuf = new void[](n + alignment);
        auto p = cast(size_t)newbuf.ptr;

        size_t i = _toNextAlignment(p) - p;
        size_t j = i + n;
        newbuf = newbuf[i .. j];

        newbuf[0 .. buf.length] = buf[];
        return newbuf;
    }
}


alias TempMemoryManager = TempMemoryManagerImpl!64;


enum bool isTempMemoryManager(T) = is(T : TempMemoryManagerImpl!alignBytes, size_t alignBytes);


unittest
{
    auto mm = new TempMemoryManagerImpl!64(128);
    static assert(isTempMemoryManager!TempMemoryManager);

    {
        auto ss1 = mm.saveState;

        auto buf1 = mm.makeArray!float(1);
        assert(((cast(size_t)buf1.ptr) & (64-1)) == 0);
        assert(mm._used == 64);

        {
            auto ss2 = mm.saveState;
            auto buf2 = mm.makeArray!float(2);
            assert(mm._used == 128);
            assert(((cast(size_t)buf2.ptr) & (64-1)) == 0);

            auto buf3 = mm.makeArray!double(1);
            assert(mm._used == 192);
            assert(((cast(size_t)buf3.ptr) & (64-1)) == 0);
        }

        assert(mm._used == 64);
    }
    assert(mm._used == 0);
}



enum hasMemoryView(T) = (is(typeof(T.init.sliced().lightScope())) && hasMemoryView!(typeof(T.init.sliced().lightScope())))
                     || is(Unqual!(typeof(T.init.sliced())) == Slice!(T.ElementType*, 1, kind), SliceKind kind)
                     || (is(Unqual!(typeof(T.init.sliced())) == Slice!(T.ElementType*, 2, kind), SliceKind kind) && (kind == Contiguous || kind == Canonical));

/** 
 * 
 */
auto makeViewOrNewSlice(Vec, Alloc)(Vec vec, ref Alloc alloc)
if(isVectorLike!Vec)
{
    static if(hasMemoryView!Vec)
    {
        return tuple!("view", "isAllocated")(vec.sliced.lightScope.lightConst, false);
    }
    else
    {
        // import std.stdio;
        // writeln("Allocate for " ~ Vec.stringof);
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


auto makeViewOrNewSlice(Mat, Alloc)(Mat mat, ref Alloc alloc)
if(isMatrixLike!Mat)
{
    static if(hasMemoryView!Mat)
    {
        return tuple!("view", "isAllocated")(mat.sliced.lightScope.lightConst, false);
    }
    else
    {
        static if(is(Unqual!(typeof(Mat.init.sliced())) == Slice!(Mat.ElementType*, 2, kind), SliceKind kind) && (kind == Universal))
        {
            auto s = mat.sliced;
            if(s._stride!0 == 1 || s._stride!1 == 1)
                return tuple!("view", "isAllocated")(s.lightConst, false);
            else
            {
                // import std.stdio;
                // writeln("Allocate Universal for " ~ Mat.stringof);
                Slice!(Mat.ElementType*, 2, Contiguous) newslice = alloc.makeArray!(Mat.ElementType)(mat.length!0 * mat.length!1).sliced(mat.length!0, mat.length!1);
                mat.evalTo(newslice, alloc);
                return tuple!("view", "isAllocated")(newslice.universal.lightConst, true);
            }
        }
        else
        {
            // import std.stdio;
            // writeln("Allocate for " ~ Mat.stringof);
            Slice!(Mat.ElementType*, 2, Contiguous) newslice = alloc.makeArray!(Mat.ElementType)(mat.length!0 * mat.length!1).sliced(mat.length!0, mat.length!1);
            mat.evalTo(newslice, alloc);
            return tuple!("view", "isAllocated")(newslice, true);
        }
    }
}

unittest
{
    auto mat1 = matrix!int(2, 2, 1);
    auto s1 = mat1.makeViewOrNewSlice(theAllocator);
    assert(s1.view == [[1, 1], [1, 1]]);
    assert(s1.isAllocated == false);

    int err;
    auto mat2 = iota(4).reshape([2, 2], err).matrixed;
    auto s2 = mat2.makeViewOrNewSlice(theAllocator);
    assert(s2.view == [[0, 1], [2, 3]]);
    assert(s2.isAllocated == true);
    theAllocator.dispose(s2.view.iterator);

    auto mat3 = mat1.T;
    auto s3 = mat3.makeViewOrNewSlice(theAllocator);
    assert(s3.view == [[1, 1], [1, 1]]);
    assert(s3.isAllocated == false);
}



struct WithMemoryView(Mat)
{
    this(Mat mat, Matrix!(Mat.ElementType, Contiguous) view)
    {
        _mat = mat;
        _view = view;
    }


    size_t length(size_t dim)() const { return _view.length!dim(); }


    Mat.ElementType opIndex(size_t i, size_t j) const
    {
        return _view[i, j];
    }


    auto sliced() const { _view.sliced.lightConst; }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        _view.evalTo(dst, alloc);
    }


  private:
    Mat _mat;
    Matrix!(Mat.ElementType, Contiguous) _view;
}


auto withMemoryView(MatA, Alloc)(Mat mat, lazy Matrix!(Mat.ElementType, Contiguous) view, ref Alloc alloc = matvecAllocator)
{
    static if(hasMemoryView!Mat)
        return mat;
    else {
        auto viewmat = view();
        mat.evalTo(viewmat.sliced, alloc);
        return WithMemoryView!Mat(mat, view);
    }
}


auto forceEvaluate(V)(V vec)
if(isVectorLike!V)
{
    auto dst = vector!(V.ElementType)(vec.length);
    dst[] = vec;
    return dst; 
}


auto forceEvaluate(M)(M mat)
if(isMatrixLike!M)
{
    auto dst = matrix!(M.ElementType)(mat.length!0, mat.length!1);
    dst[] = mat;
    return dst; 
}


unittest
{
    Vector!(float, Contiguous) vec = vector!float(3, 0);
    vec[0] = 1;
    vec[1] = 2;
    vec[2] = 3;

    Vector!(float, Contiguous) cache = forceEvaluate(vec * 2 + vec);
    assert(cache[0] == 3);
    assert(cache[1] == 6);
    assert(cache[2] == 9);
}


// 
struct MatrixMatrixMulGEMM(S, MatA, MatB, MatC)
{
    alias ElementType = typeof(S.init * MatA.ElementType.init * MatB.ElementType.init + S.init * MatC.ElementType.init);
    enum float exprTreeCost = sum([staticMap!(.exprTreeCost, AliasSeq!(MatA, MatB, MatC))]) + EXPR_COST_N^^3;


    this(S alpha, MatA matA, MatB matB, S beta, MatC matC)
    {
        _alpha = alpha;
        _beta = beta;
        _matA = matA;
        _matB = matB;
        _matC = matC;
    }


    size_t length(size_t dim)() const
    if(dim == 0 || dim == 1)
    {
        static if(dim == 0)
            return _matA.length!0;
        else
            return _matB.length!1;
    }


    ElementType opIndex(size_t i, size_t j) const
    {
        ElementType sum = ElementType(0);
        foreach(n; 0 .. _matA.length!1)
            sum += _matA[i, n] * _matB[n, j];
        
        return sum * _alpha + _beta * _matC[i, j];
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc)
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        static if(MatA.exprTreeCost < EXPR_COST_N && MatB.exprTreeCost < EXPR_COST_N && MatC.exprTreeCost < EXPR_COST_N)
        {
            if(_matA.length!0 * _matA.length!1 * _matB.length!1 <= 100) {
                _evalToForSmall(dst);
                return;
            }
        }

        auto viewA = _matA.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        auto viewB = _matB.makeViewOrNewSlice(alloc);
        scope(exit) if(viewB.isAllocated) alloc.dispose(viewB.view.iterator);

        _matC.evalTo(dst, alloc);

        static if(isBlasType!T)
        {
            static if(isFloatingPoint!T || (is(T == MirComplex!E, E) && isFloatingPoint!E))
            {
                import mir.blas : gemm;
                gemm(_alpha, viewA.view, viewB.view, _beta, dst);
            }
            else
            {
                import dffdd.math.linalg : gemm_stdcomplex;
                gemm_stdcomplex(T(_alpha), viewA.view, viewB.view, T(_beta), dst);
            }
        }
        else
        {
            import dffdd.math.linalg : dot;
            auto trB = viewB.view.lightConst.transposed;
            foreach(i; 0 .. this.length!0)
                foreach(j; 0 .. this.length!1) {
                    dst[i, j] *= _beta;
                    dst[i, j] += dot(viewA.view[i].vectored, trB[j].vectored) * _alpha;
                }
        }
    }


    mixin(definitionsOfMatrixOperators(
            ["defaults", "M*M"]
            ~ (isZeros!MatC ? [] : ["M+M"])));


    auto T()()
    {
        return _alpha * (_matB.T * _matA.T) + _beta * _matC.T;
    }


  static if(isComplex!ElementType)
  {
    auto H()()
    {
        return conj(_alpha) * (_matB.H * _matA.H) + conj(_beta) * _matC.H;
    }
  }
  else
  {
    alias H = T;
  }


    auto _opB_Impl_(string op : "*", U)(U u)
    if(!isMatrixLike!U && !isVectorLike!U && is(typeof(U.init * S.init)))
    {
        return matrixGemm(_alpha * u, _matA, _matB, _beta * u, _matC);
    }


    auto _opBR_Impl_(string op : "*", U)(U u)
    if(!isMatrixLike!U && !isVectorLike!U && is(typeof(U.init * S.init)))
    {
        return matrixGemm(_alpha * u, _matA, _matB, _beta * u, _matC);
    }


    auto _opB_Impl_(string op : "*", V)(V vec)
    if(isVectorLike!V)
    in(vec.length == this.length!1)
    {
        static if(isZeros!MatC)
        {
            return _matA * (_matB * vec) * _alpha;
        }
        else
        {
            auto vec2 = forceEvaluate(vec);
            return _matA * (_matB * vec2) * _alpha + (_matC * vec2) * _beta;
        }
    }


  static if(isZeros!MatC)
  {
    auto _opB_Impl_(string op, M)(M mat)
    if(isMatrixLike!M && (op == "+" || op == "-"))
    in(mat.length!0 == this.length!0 && mat.length!1 == this.length!1)
    {
        static if(is(M == MatrixAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
            return matrixGemm(_alpha, _matA, _matB, op == "+" ? mat._alpha : -mat._alpha, mat._matA);
        else
            return matrixGemm(_alpha, _matA, _matB, ElementType(op == "+" ? 1 : -1), mat);
    }


    auto _opBR_Impl_(string op, M)(M mat)
    if(isMatrixLike!M && (op == "+" || op == "-"))
    in(mat.length!0 == this.length!0 && mat.length!1 == this.length!1)
    {
        static if(is(M == MatrixAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
            return matrixGemm(op == "+" ? _alpha : -_alpha, _matA, _matB, mat._alpha, mat._matA);
        else
            return matrixGemm(op == "+" ? _alpha : -alpha, _matA, _matB, ElementType(1), mat);
    }
  }
//   else
//   {
//     mixin(definitionsOfMatrixOperators(["M+M"]));
//   }


  private:
    S _alpha;
    S _beta;
    MatA _matA;
    MatB _matB;
    MatC _matC;


    void _evalToForSmall(T, SliceKind kindA)(Slice!(T*, 2, kindA) dst) const
    {
        gemm_smallsize(_alpha, _matA, _matB, _beta, _matC, dst);
    }
}


private
void gemm_smallsize(T, MatA, MatB, MatC, SliceKind kindD)
    (T alpha, const ref MatA a, const ref MatB b, T beta, const ref MatC c, ref Slice!(T*, 2, kindD) dst)
in
{
    assert(a.length!1 == b.length!0);
    assert(a.length!0 == c.length!0);
    assert(c.length!1 == b.length!1);
    assert(c.length!0 == dst.length!0);
}
do
{
    immutable M = a.length!0,
              N = b.length!1,
              K = a.length!1;

    dst[] = 0;

    static if(isZeros!MatC)
    {
        foreach(i; 0 .. M) {
            foreach(k; 0 .. K)
                foreach(j; 0 .. N)
                    dst[i, j] += a[i, k] * b[k, j];

            foreach(j; 0 .. N)
                dst[i, j] = dst[i, j] * alpha + beta * c[i, j];
        }
    }
    else
    {
        foreach(i; 0 .. M)
            foreach(j; 0 .. N) {
                T s = T(0);
                foreach(k; 0 .. K)
                    s += a[i, k] * b[k, j];

                dst[i, j] = s * alpha + beta * c[i, j];
            }
    }
}


MatrixMatrixMulGEMM!(T, MatA, MatB, MatC) matrixGemm(T, MatA, MatB, MatC)(T alpha, MatA matA, MatB matB, T beta, MatC matC)
{
    return typeof(return)(alpha, matA, matB, beta, matC);
}


unittest
{
    static foreach(E; AliasSeq!(int, long, float, double, StdComplex!float, MirComplex!double))
    {{
        int err;
        auto mat1 = iota(6).reshape([2, 3], err).as!E.matrixed;
        auto mat2 = iota(6).reshape([3, 2], err).as!E.matrixed;
        auto mat3 = iota(4).reshape([2, 2], err).as!E.matrixed;

        auto mul1 = matrixGemm(E(2), mat1, mat2, E(3), mat3);
        static assert(isMatrixLike!(typeof(mul1)));
        assert(err == 0);

        assert(mul1.length!0 == 2);
        assert(mul1.length!1 == 2);
        assert(mul1[0, 0] == 20);
        assert(mul1[0, 1] == 29);
        assert(mul1[1, 0] == 62);
        assert(mul1[1, 1] == 89);

        auto s = slice!E(2, 2);
        mul1.evalTo(s, theAllocator);
        assert(s[0, 0] == 20);
        assert(s[0, 1] == 29);
        assert(s[1, 0] == 62);
        assert(s[1, 1] == 89);

        auto view = mul1.makeViewOrNewSlice(theAllocator);
        assert(view.view == [[20, 29], [62, 89]]);
        assert(view.isAllocated == true);
        static assert(isMatrixLike!(typeof(mul1)));
    }}
}

unittest
{
    auto mat1 = [0, 1, 2, 3, 4, 5].sliced(2, 3).matrixed;
    auto mat2 = [0, 1, 2, 3, 4, 5].sliced(3, 2).matrixed;
    auto mat3 = [0, 1, 2, 3].sliced(2, 2).matrixed;

    auto vec1 = [2, 3].sliced.vectored;

    // auto z1 = mat1 * mat2 * 2 + mat3 * 3;
    auto z1 = (mat1 * mat2 * 2)._opB_Impl_!"+"(mat3 * 3);
    assert(z1.makeViewOrNewSlice(theAllocator).view == [[20, 29], [62, 89]]);

    auto z2 = z1 * vec1;
    assert(z2.makeViewOrNewSlice(theAllocator).view == [127, 391]);
}




// 
struct MatrixVectorMulGEMV(T, MatA, VecB, VecC)
if(isMatrixLike!MatA && isVectorLike!VecB && !isMatrixLike!T && !isVectorLike!T)
{
    alias ElementType = typeof(T.init * MatA.ElementType.init * VecB.ElementType.init + T.init * VecC.ElementType.init);
    enum float exprTreeCost = sum([staticMap!(.exprTreeCost, AliasSeq!(MatA, VecB, VecC))]) + EXPR_COST_N^^2;


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


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc)
    in(dst.length == this.length)
    {
        auto viewA = _matA.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        auto viewB = _vecB.makeViewOrNewSlice(alloc);
        scope(exit) if(viewB.isAllocated) alloc.dispose(viewB.view.iterator);

        _vecC.evalTo(dst, alloc);

        static if(isBlasType!T && _useGEMV)
        {
            static if(isFloatingPoint!T || (is(T == MirComplex!E, E) && isFloatingPoint!E))
            {
                import mir.blas : gemv;
                gemv(_alpha, viewA.view, viewB.view, _beta, dst);
            }
            else
            {
                import dffdd.math.linalg : gemv_stdcomplex;
                gemv_stdcomplex(T(_alpha), viewA.view, viewB.view, T(_beta), dst);
            }
        }
        else
        {
            import dffdd.math.linalg : dot;            

            foreach(n; 0 .. this.length) {
                dst[n] *= _beta;
                dst[n] += dot(viewA.view[n].vectored, viewB.view.vectored) * _alpha;
            }
        }
    }


    mixin(definitionsOfVectorOperators(["defaults"]));


    auto _opB_Impl_(string op : "*", S)(S val)
    if(is(T :  ElementType))
    {
        return vectorGemv(_alpha * val, _matA, _vecB, _beta * val, _vecC);
    }


    auto _opBR_Impl_(string op : "*", S)(S val)
    if(is(T : ElementType))
    {
        return vectorGemv(_alpha * val, _matA, _vecB, _beta * val, _vecC);
    }


    static if(isZeros!VecC)
    {
        auto _opB_Impl_(string op, V)(V vec)
        if(isVectorLike!V && (op == "+" || op == "-"))
        in(vec.length == this.length)
        {
            static if(is(V == VectorAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
                return vectorGemv(_alpha, _matA, _vecB, op == "+" ? vec._alpha : -vec._alpha, vec._vecA);
            else
                return vectorGemv(_alpha, _matA, _vecB, ElementType(op == "+" ? 1 : -1), vec);
        }

        auto _opBR_Impl_(string op, V)(V vec)
        if(isVectorLike!V && (op == "+" || op == "-"))
        in(vec.length == this.length)
        {
            static if(is(V == VectorAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
                return vectorGemv(op == "+" ? _alpha : -_alpha, _matA, _vecB, vec._alpha, vec._vecA);
            else
                return vectorGemv(op == "+" ? _alpha : -_alpha, _matA, _vecB, ElementType(1), vec);
        }
    }
    else
    {
        mixin(definitionsOfVectorOperators(["V+V"]));
    }


  private:
    T _alpha;
    T _beta;
    MatA _matA;
    VecB _vecB;
    VecC _vecC;


    version(DFFDD_BLAS_USEGEMV)
        enum bool _useGEMV = true;
    else
        enum bool _useGEMV = false;
}


MatrixVectorMulGEMV!(T, MatA, VecB, VecC) vectorGemv(T, MatA, VecB, VecC)(T alpha, MatA matA, VecB vecB, T beta, VecC vecC)
{
    return typeof(return)(alpha, matA, vecB, beta, vecC);
}

unittest
{
    static foreach(E; AliasSeq!(int, long, float, double, StdComplex!float, MirComplex!double))
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


// aX + bY
struct MatrixAxpby(S, MatA, MatB)
if(isMatrixLike!MatA && isMatrixLike!MatB)
{
    alias ElementType = typeof(S.init * MatA.ElementType.init + S.init * MatB.ElementType.init);
    enum float exprTreeCost = sum([staticMap!(.exprTreeCost, AliasSeq!(MatA, MatB))]) + EXPR_COST_N^^2;


    this(S alpha, MatA matA, S beta, MatB matB)
    {
        _alpha = alpha;
        _beta = beta;
        _matA = matA;
        _matB = matB;
    }


    size_t length(size_t dim)() const
    {
        return _matA.length!dim;
    }


    ElementType opIndex(size_t i, size_t j) const
    {
        return _alpha * _matA[i, j] + _beta * _matB[i, j];
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc)
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        static if(!hasMemoryView!MatA && hasMemoryView!MatB)
        {
            matrixAxpby(_beta, _matB, _alpha, _matA).evalTo(dst, alloc);
        }
        else
        {
            auto viewA = _matA.makeViewOrNewSlice(alloc);
            scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

            _matB.evalTo(dst, alloc);
            dst[] *= _beta;

            foreach(i; 0 .. this.length!0)
                foreach(j; 0 .. this.length!1)
                    dst[i, j] += _alpha * viewA.view[i, j];
        }
    }


    // mixin MatrixOperators!(["M+M", "M*M"]) OpImpls;
    mixin(definitionsOfMatrixOperators(
        ["defaults", ".H", ".T"]
        ~ (isZeros!MatB ? [] : ["M+M", "M*M", "M*V"])
    ));


    auto _opB_Impl_(string op : "*", U)(U u)
    if(!isVectorLike!U && !isMatrixLike!U && is(typeof(U.init * S.init)))
    {
        return matrixAxpby(_alpha * u, _matA, _beta * u, _matB);
    }


    auto _opBR_Impl_(string op : "*", U)(U u)
    if(!isVectorLike!U && !isMatrixLike!U && is(typeof(U.init * S.init)))
    {
        return matrixAxpby(_alpha * u, _matA, _beta * u, _matB);
    }


  static if(isZeros!MatB)
  {
    auto _opB_Impl_(string op, M)(M mat)
    if(isMatrixLike!M && (op == "+" ||  op == "-"))
    in(mat.length!0 == this.length!0 && mat.length!1 == this.length!1)
    {
        static if(is(M == MatrixAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
            return matrixAxpby(_alpha, _matA, op == "+" ? mat._alpha : -mat._alpha, mat._matA);
        else
            return matrixAxpby(_alpha, _matA, S(op == "+" ? 1 : -1), mat);
    }


    auto _opB_Impl_(string op : "*", V)(V vec)
    if(isVectorLike!V)
    in(vec.length == this.length!1)
    {
        return _alpha * (_matA * vec);
    }


    auto _opB_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(mat.length!0 == this.length!1)
    {
        return _alpha * (_matA * mat);
    }


    auto _opBR_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(mat.length!1 == this.length!0)
    {
        return _alpha * (mat * _matA);
    }
  }
//   else
//   {
//     mixin(definitionsOfMatrixOperators(["M+M", "M*M", "M*V"]));
//   }


  private:
    S _alpha;
    S _beta;
    MatA _matA;
    MatB _matB;
}


MatrixAxpby!(T, MatA, MatB) matrixAxpby(T, MatA, MatB)(T scalarA, MatA vecA, T scalarB, MatB vecB)
{
    return typeof(return)(scalarA, vecA, scalarB, vecB);
}

unittest
{
    static foreach(E; AliasSeq!(int, long, float, double, StdComplex!float, MirComplex!double))
    {{
        int err;
        auto mat1 = iota(6).reshape([2, 3], err).as!E.matrixed;
        auto mat2 = 4.repeat([2, 3]).as!E.matrixed;

        auto add1 = matrixAxpby(E(2), mat1, E(1), mat2);
        static assert(isMatrixLike!(typeof(add1)));
        assert(err == 0);

        assert(add1.length!0 == 2);
        assert(add1.length!1 == 3);
        assert(add1[0, 0] == 4);
        assert(add1[0, 1] == 6);
        assert(add1[0, 2] == 8);
        assert(add1[1, 0] == 10);
        assert(add1[1, 1] == 12);
        assert(add1[1, 2] == 14);

        auto s = slice!E(2, 3);
        add1.evalTo(s, theAllocator);
        assert(s[0, 0] == 4);
        assert(s[0, 1] == 6);
        assert(s[0, 2] == 8);
        assert(s[1, 0] == 10);
        assert(s[1, 1] == 12);
        assert(s[1, 2] == 14);
    }}
}



// aX + bY
struct VectorAxpby(T, VecA, VecB)
if(isVectorLike!VecA && isVectorLike!VecB)
{
    alias ElementType = typeof(T.init * VecA.ElementType.init + T.init * VecB.ElementType.init);
    enum float exprTreeCost = sum([staticMap!(.exprTreeCost, AliasSeq!(VecA, VecB))]) + EXPR_COST_N;


    this(T alpha, VecA vecA, T beta, VecB vecB)
    {
        _alpha = alpha;
        _beta = beta;
        _vecA = vecA;
        _vecB = vecB;
    }


    ElementType opIndex(size_t n) const
    in(n < this.length)
    {
        return _alpha * _vecA[n] + _beta * _vecB[n];
    }


    size_t length() const
    {
        return _vecA.length;
    }


    void evalTo(E, SliceKind kindA, Alloc)(Slice!(E*, 1, kindA) dst, ref Alloc alloc)
    in(dst.length == this.length)
    {
        static if(!hasMemoryView!VecA && hasMemoryView!VecB)
        {
            vectorAxpby(_beta, _vecB, _alpha, _vecA).evalTo(dst, alloc);
        }
        else
        {
            auto viewA = _vecA.makeViewOrNewSlice(alloc);
            scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

            _vecB.evalTo(dst, alloc);
            dst[] *= _beta;

            static if(isFloatingPoint!E || (is(E == MirComplex!F, F) && isFloatingPoint!F && is(E == VecA.ElementType) ))
            {
                import mir.blas : axpy;
                static if(is(E == T))
                    axpy(_alpha, viewA.view, dst);
                else
                    axpy(E(_alpha), viewA.view, dst);
            }
            else
            {
                foreach(i; 0 .. this.length)
                    dst[i] += _alpha * viewA.view[i];
            }
        }
    }


    mixin VectorOperators!(["defaults", "V+V"]) OpImpls;


    auto _opB_Impl_(string op : "*", U)(U u)
    if(!isVectorLike!U && !isMatrixLike!U && is(typeof(U.init * T.init)))
    {
        return vectorAxpby(_alpha * u, _vecA, _beta * u, _vecB);
    }


    auto _opBR_Impl_(string op : "*", U)(U u)
    if(!isVectorLike!U && !isMatrixLike!U && is(typeof(U.init * T.init)))
    {
        return vectorAxpby(_alpha * u, _vecA, _beta * u, _vecB);
    }


    auto _opB_Impl_(string op : "/", U)(U u)
    if(!isVectorLike!U && !isMatrixLike!U && is(typeof(T.init / U.init)))
    {
        return vectorAxpby(_alpha / u, _vecA, _beta / u, _vecB);
    }


    static if(isZeros!VecB)
    {
        auto _opB_Impl_(string op, V)(V vec)
        if(isVectorLike!V && (op == "+" || op == "-"))
        in(vec.length == vec.length)
        {
            static if(is(V == VectorAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
                return vectorAxpby(_alpha, _vecA, op == "+" ? vec._alpha : -vec._alpha, vec._vecA);
            else
                return vectorAxpby(_alpha, _vecA, T(op == "+" ? 1 : -1), vec);
        }

        auto _opBR_Impl_(string op, V)(V vec)
        if(isVectorLike!V && (op == "+" || op == "-"))
        in(vec.length == vec.length)
        {
            static if(is(V == VectorAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
                return vectorAxpby(op == "+" ? _alpha : -_alpha, _vecA, vec._alpha, vec._vecA);
            else
                return vectorAxpby(op == "+" ? _alpha : -_alpha, _vecA, T(1), vec);
        }
    }
    else
    {
        mixin(definitionsOfVectorOperators(["V+V"]));
    }


  private:
    T _alpha;
    T _beta;
    VecA _vecA;
    VecB _vecB;
}


VectorAxpby!(T, VecA, VecB) vectorAxpby(T, VecA, VecB)(T scalarA, VecA vecA, T scalarB, VecB vecB)
{
    return typeof(return)(scalarA, vecA, scalarB, vecB);
}


unittest
{
    static foreach(E; AliasSeq!(int, long, float, double, StdComplex!float, MirComplex!double))
    {{
        E alpha = E(2);
        auto vec1 = iota(3).as!E.vectored;
        auto vec2 = [3, 2, 3].sliced.as!E.vectored;

        auto axpy1 = vectorAxpby(alpha, vec1, E(1), vec2);
        axpy1.evalTo(slice!E(3), theAllocator);
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


struct ElementMap(alias fn, Mat)
if(isMatrixLike!Mat || isVectorLike!Mat)
{
    import std.functional : unaryFun;
    alias ElementType = typeof(unaryFun!fn(Mat.ElementType.init));

  static if(isMatrixLike!Mat)
  {
    enum float exprTreeCost = Mat.exprTreeCost + EXPR_COST_N^^2;
  }
  else
  {
    enum float exprTreeCost = Mat.exprTreeCost + EXPR_COST_N;
  }


    this(Mat mat)
    {
        _mat = mat;
    }


  static if(isMatrixLike!Mat)
  {
    size_t length(size_t dim)() const { return _mat.length!dim(); }


    ElementType opIndex(size_t i, size_t j) const
    {
        return unaryFun!fn(_mat[i, j]);
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc)
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        auto viewA = _mat.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        dst[] = viewA.view.map!fn;
    }


    mixin MatrixOperators!(["defaults", "M*M", "M*V", "M+M", ".T", ".H"]);
  }
  else
  {
    size_t length() const { return _mat.length; }


    ElementType opIndex(size_t i) const
    {
        return fn(_mat[i, j]);
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc)
    in(dst.length == this.length)
    {
        auto viewA = _mat.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        dst[] = _mat.view.map!fn;
    }


    mixin VectorOperators!(["defaults", "S*V", "V*S", "V+V"]);
  }


  private:
    Mat _mat;
}

auto elementMap(alias fn, Mat)(Mat mat)
if(isMatrixLike!Mat || isVectorLike!Mat)
{
    return ElementMap!(fn, Mat)(mat);
}

unittest
{
    int err;
    auto mat = iota(4).reshape([2, 2], err).matrixed.elementMap!"a*2";
    // static assert(isMatrixLike!(typeof(mat)));

    assert(mat[0, 0] == 0);
    assert(mat[0, 1] == 2);
    assert(mat[1, 0] == 4);
    assert(mat[1, 1] == 6);

    auto fn = (typeof(mat) t, size_t i){
        alias T = typeof(t);
        T.ElementType e = t[i, i];
        size_t rowlen = t.length!0;
        size_t collen = t.length!1;

        auto dst = slice!(typeof(e))(rowlen, collen);
        auto allocator = std.experimental.allocator.theAllocator;
        t.evalTo(dst, allocator);
    };
}


struct TransposedAndOrConjugated(Flag!"isTransposed" isTransposed, Flag!"isConjugated" isConjugated, Mat)
if(isConjugated ? isComplex!(Mat.ElementType) : true)
{
    alias ElementType = Mat.ElementType;
    enum float exprTreeCost = Mat.exprTreeCost + 1;


    this(Mat mat)
    {
        _mat = mat;
    }


    ElementType opIndex(size_t i, size_t j) const
    {
        ElementType ret;
        static if(isTransposed)
            ret = _mat[j, i];
        else
            ret = _mat[i, j];
        
        static if(isConjugated && isComplex!ElementType)
            return conj(ret);
        else
            return ret;
    }


    size_t length(size_t dim)() const
    if(dim == 0 || dim == 1)
    {
        static if((dim == 0 && isTransposed) || (dim == 1 && !isTransposed))
            return _mat.length!1;
        else
            return _mat.length!0;
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc)
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        static if(isTransposed)
            _mat.evalTo(dst.transposed, alloc);
        else
            _mat.evalTo(dst, alloc);

        static if(isConjugated && isComplex!ElementType) {
            foreach(i; 0 .. dst.length!0)
                foreach(j; 0 .. dst.length!1)
                    dst[i, j] = conj(dst[i, j]);
        }
    }


    auto T()
    {
        static if(isTransposed && !isConjugated)
            return _mat;
        else
            return TransposedAndOrConjugated!(isTransposed ? No.isTransposed : Yes.isTransposed, isConjugated, Mat)(_mat);
    }


  static if(isComplex!ElementType)
  {
    auto H()
    {
        static if(isTransposed && isConjugated)
            return _mat;
        else
            return TransposedAndOrConjugated!(isTransposed ? No.isTransposed : Yes.isTransposed, isConjugated ? No.isConjugated : Yes.isConjugated, Mat)(_mat);
    }
  }
  else
  {
    alias H = T;
  }


    mixin MatrixOperators!(["defaults", "M*M", "M*V", "M+M"]);


  private:
    Mat _mat;
}


TransposedAndOrConjugated!(Yes.isTransposed, No.isConjugated, M) matrixTransposed(M)(M mat)
{
    return typeof(return)(mat);
}


unittest
{
    static foreach(E; AliasSeq!(int, long, float, double, StdComplex!float, MirComplex!double))
    {{
        int err;
        auto mat1 = iota(6).reshape([3, 2], err).as!E.matrixed;
        auto matT = matrixTransposed(mat1);
        assert(matT.length!0 == 2);
        assert(matT.length!1 == 3);
        static assert(isMatrixLike!(typeof(mat1)));

        assert(matT[0, 0] == 0);
        assert(matT[0, 1] == 2);
        assert(matT[0, 2] == 4);
        assert(matT[1, 0] == 1);
        assert(matT[1, 1] == 3);
        assert(matT[1, 2] == 5);
    }}
}


mixin template VectorOperators(string[] list)
{
    mixin(definitionsOfVectorOperators(list));
}


string definitionsOfVectorOperators(string[] list)
{
    import std.algorithm;

    string dst = q{
    import std.traits : isFloatingPoint, isIntegral;
    import dffdd.math.linalg : isBlasType;
    import dffdd.math.exprtemplate;
    import dffdd.math.matrix;
    import dffdd.math.vector;
    import dffdd.math.matrixspecial;
    };


    if(canFind(list, "defaults"))
    dst ~= q{
        unittest {
            typeof(this) a;
            checkIsVectorLike(a);
        }

        auto opBinary(string op, X)(X rhs)
        if(op == "+" || op == "-" || op == "*")
        {
            static if(is(typeof(this._opB_Impl_!op(rhs))) && is(typeof(rhs._opBR_Impl_!op(this))))
            {
                alias T1 = typeof(this._opB_Impl_!op(rhs));
                alias T2 = typeof(rhs._opBR_Impl_!op(this));
                static if(T1.exprTreeCost <= T2.exprTreeCost)
                    return this._opB_Impl_!op(rhs);
                else
                    return rhs._opBR_Impl_!op(this);
            }
            else static if(is(typeof(this._opB_Impl_!op(rhs))))
            {
                return this._opB_Impl_!op(rhs);
            }
            else static if(is(typeof(rhs._opBR_Impl_!op(this))))
            {
                return rhs._opBR_Impl_!op(this);
            }
            else static assert(0, "Operator '(LHS) " ~ op ~ " (RHS)' is not defined with (LHS) = " ~ typeof(this).stringof ~ ", and (RHS) = " ~ X.stringof);
        }


        auto opBinary(string op : "/", X)(X rhs)
        if(!isMatrixLike!X && !isVectorLike!X && is(typeof(ElementType.init / rhs)))
        {
            static if(is(typeof(this._opB_Impl_!"/"(rhs))))
                return this._opB_Impl_!"/"(rhs);
            else
                return this.opBinary!"*"(rhs^^(-1));
        }


        auto opBinaryRight(string op, X)(X rhs)
        if(op == "+" || op == "-" || op == "*")
        {
            static if(is(typeof(this._opBR_Impl_!op(rhs))) && is(typeof(rhs._opB_Impl_!op(this))))
            {
                alias T1 = typeof(this._opBR_Impl_!op(rhs));
                alias T2 = typeof(rhs._opB_Impl_!op(this));
                static if(T1.exprTreeCost <= T2.exprTreeCost)
                    return this._opBR_Impl_!op(rhs);
                else
                    return rhs._opB_Impl_!op(this);
            }
            else static if(is(typeof(this._opBR_Impl_!op(rhs))))
            {
                return this._opBR_Impl_!op(rhs);
            }
            else static if(is(typeof(rhs._opB_Impl_!op(this))))
            {
                return rhs._opB_Impl_!op(this);
            }
            else static if(!isMatrixLike!X && !isVectorLike!X && is(X : ElementType) && (op == "*" || op == "+"))
            {
                return this.opBinary!op(rhs);
            }
            else static if(!isMatrixLike!X && !isVectorLike!X && is(X : ElementType) && (op == "-"))
            {
                return this.opBinary!"-"(-rhs) * ElementType(-1);
            }
            else static assert(0, "Operator '(LHS) " ~ op ~ " (RHS)' is not defined with (LHS) = " ~ X.stringof ~ ", and (RHS) = " ~ typeof(this).stringof);
        }
    };


    if(canFind(list, "V*S") || canFind(list, "S*V"))
    dst ~= q{
        auto _opB_Impl_(string op : "*", S)(S scalar)
        if(!isVectorLike!S && !isMatrixLike!S && is(typeof(S.init * ElementType.init)))
        {
            return vectorAxpby(scalar, this, S(0), zeros!ElementType(this.length));
        }

        auto _opBR_Impl_(string op : "*", S)(S scalar)
        if(!isVectorLike!S && !isMatrixLike!S && is(typeof(S.init * ElementType.init)))
        {
            return vectorAxpby(scalar, this, S(0), zeros!ElementType(this.length));
        }

        auto _opB_Impl_(string op : "/", S)(S scalar)
        if(!isVectorLike!S && !isMatrixLike!S && is(typeof(ElementType.init / S.init)))
        {
            return vectorAxpby(1/scalar, this, S(0), zeros!ElementType(this.length));
        }
    };


    if(canFind(list, "V+V"))
    dst ~= q{
        auto _opB_Impl_(string op, V)(V vec)
        if(isVectorLike!V && (op == "+" || op == "-"))
        in(vec.length == this.length)
        {
            static if(isZeros!V)
                return this;
            else
                return vectorAxpby(ElementType(1), this, ElementType(op == "+" ? 1 : -1), vec);
        }
    };


    if(canFind(list, "V=V"))
    dst ~= q{
        auto opSliceAssign(V)(V vec)
        if(isVectorLike!V)
        in(vec.length == this.length)
        {
            auto view = vec.makeViewOrNewSlice(matvecAllocator);
            scope(exit) if(view.isAllocated) matvecAllocator.dispose(view.view.iterator);

            foreach(n; 0 .. this.length)
                this[n] = view.view[n];
        }
    };

    if(canFind(list, "V=S"))
    dst ~= q{
        auto opSliceAssign(T)(T scalar)
        if(is(typeof(ElementType(T.init))))
        {
            foreach(n; 0 .. this.length)
                this[n] = scalar;
        }
    };


    return dst;
}


mixin template MatrixOperators(string[] list)
{
    mixin(definitionsOfMatrixOperators(list));
}


string definitionsOfMatrixOperators(string[] list)
{
    import std.algorithm : canFind;

    string dst = q{
    import dffdd.math.complex : isComplex;
    import std.typecons : Tuple, Flag, Yes, No;
    import std.traits : isFloatingPoint, isIntegral;
    import dffdd.math.linalg : isBlasType;
    import dffdd.math.exprtemplate;
    import dffdd.math.matrix;
    import dffdd.math.vector;
    import dffdd.math.matrixspecial;
    };


    if(canFind(list, "defaults"))
    dst ~= q{
        unittest {
            typeof(this) a;
            checkIsMatrixLike(a);
        }


        auto opBinary(string op, X)(X rhs)
        if(op == "+" || op == "-" || op == "*")
        {
            static if(is(typeof(this._opB_Impl_!op(rhs))) && is(typeof(rhs._opBR_Impl_!op(this))))
            {
                alias T1 = typeof(this._opB_Impl_!op(rhs));
                alias T2 = typeof(rhs._opBR_Impl_!op(this));
                static if(T1.exprTreeCost <= T2.exprTreeCost)
                    return this._opB_Impl_!op(rhs);
                else
                    return rhs._opBR_Impl_!op(this);
            }
            else static if(is(typeof(this._opB_Impl_!op(rhs))))
            {
                return this._opB_Impl_!op(rhs);
            }
            else static if(is(typeof(rhs._opBR_Impl_!op(this))))
            {
                return rhs._opBR_Impl_!op(this);
            }
            else static assert(0, "Operator '(LHS) " ~ op ~ " (RHS)' is not defined with (LHS) = " ~ typeof(this).stringof ~ ", and (RHS) = " ~ X.stringof);
        }


        auto opBinary(string op : "/", X)(X rhs)
        if(!isMatrixLike!X && !isVectorLike!X && is(typeof(ElementType.init / rhs)))
        {
            static if(is(typeof(this._opB_Impl_!"/"(rhs))))
                return this._opB_Impl_!"/"(rhs);
            else
                return this.opBinary!"*"(rhs^^(-1));
        }


        auto opBinaryRight(string op, X)(X rhs)
        if(op == "+" || op == "-" || op == "*")
        {
            static if(is(typeof(this._opBR_Impl_!op(rhs))) && is(typeof(rhs._opB_Impl_!op(this))))
            {
                alias T1 = typeof(this._opBR_Impl_!op(rhs));
                alias T2 = typeof(rhs._opB_Impl_!op(this));
                static if(T1.exprTreeCost <= T2.exprTreeCost)
                    return this._opBR_Impl_!op(rhs);
                else
                    return rhs._opB_Impl_!op(this);
            }
            else static if(is(typeof(this._opBR_Impl_!op(rhs))))
            {
                return this._opBR_Impl_!op(rhs);
            }
            else static if(is(typeof(rhs._opB_Impl_!op(this))))
            {
                return rhs._opB_Impl_!op(this);
            }
            else static if(!isMatrixLike!X && !isVectorLike!X && is(X : ElementType) && (op == "*" || op == "+"))
            {
                return this.opBinary!op(rhs);
            }
            else static if(!isMatrixLike!X && !isVectorLike!X && is(X : ElementType) && (op == "-"))
            {
                return this.opBinary!"-"(-rhs) * ElementType(-1);
            }
            else static assert(0, "Operator '(LHS) " ~ op ~ " (RHS)' is not defined with (LHS) = " ~ X.stringof ~ ", and (RHS) = " ~ typeof(this).stringof);
        }


        size_t opDollar(size_t dim)()
        {
            return this.length!dim;
        }


        Tuple!(size_t, size_t) opSlice(size_t dim)(size_t i, size_t j)
        {
            return typeof(return)(i, j);
        }


        auto opIndex(Tuple!(size_t, size_t) rs, Tuple!(size_t, size_t) cs)
        {
            static if(is(typeof(this._op_partial_impl_(rs, cs))))
            {
                return this._op_partial_impl_(rs, cs);
            }
            else
            {
                return PartialFromMatrix!(typeof(this), typeof(rs), typeof(cs))(this, rs, cs);
            }
        }


        auto opIndex(Tuple!(size_t, size_t) rs, size_t j)
        {
            static if(is(typeof(this._op_partial_impl_(rs, j))))
            {
                return this._op_partial_impl_(rs, j);
            }
            else
            {
                return PartialFromMatrix!(typeof(this), typeof(rs), size_t)(this, rs, j);
            }
        }


        auto opIndex(size_t i, Tuple!(size_t, size_t) cs)
        {
            static if(is(typeof(this._op_partial_impl_(i, cs))))
            {
                return this._op_partial_impl_(i, cs);
            }
            else
            {
                return PartialFromMatrix!(typeof(this), size_t, typeof(cs))(this, i, cs);
            }
        }
    };


    if(canFind(list, "M*M"))
    dst ~= q{
        auto _opB_Impl_(string op : "*", M)(M mat)
        if(isMatrixLike!M)
        in(mat.length!0 == this.length!1)
        {
            static if(isIdentity!M)
                return this;
            else static if(isZeros!M)
                return zeros!ElementType(this.length!0, mat.length!1);
            else
                return matrixGemm(ElementType(1), this, mat, ElementType(0), zeros!ElementType(this.length!0, mat.length!1));
        }
    };


    if(canFind(list, "M*S") || canFind(list, "S*M"))
    dst ~= q{
        auto _opB_Impl_(string op : "*", S)(S scalar)
        if(!isVectorLike!S && !isMatrixLike!S && is(typeof(S.init * ElementType.init)))
        {
            return matrixAxpby(scalar, this, S(0), zeros!ElementType(this.length!0, this.length!1));
        }


        auto _opBR_Impl_(string op : "*", S)(S scalar)
        if(!isVectorLike!S && !isMatrixLike!S && is(typeof(S.init * ElementType.init)))
        {
            return matrixAxpby(scalar, this, S(0), zeros!ElementType(this.length!0, this.length!1));
        }
    };


    if(canFind(list, "M+M"))
    dst ~= q{
        auto _opB_Impl_(string op, M)(M mat)
        if(isMatrixLike!M && (op == "+" || op =="-"))
        in(mat.length!0 == this.length!0 && mat.length!1 == this.length!1)
        {
            static if(isZeros!M)
                return this;
            else {
                static if(is(M == MatrixAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
                    return matrixAxpby(ElementType(1), this, op == "+" ? mat._alpha : -mat._alpha, mat._matA);
                else
                    return matrixAxpby(ElementType(1), this, ElementType(op == "+" ? 1 : -1), mat);
            }
        }
    };


    if(canFind(list, "M*V"))
    dst ~= q{
        auto _opB_Impl_(string op : "*", V)(V vec)
        if(isVectorLike!V)
        in(vec.length == this.length!1)
        {
            return vectorGemv(ElementType(1), this, vec, ElementType(0), zeros!ElementType(this.length!0));
        }
    };

    if(canFind(list, ".T"))
    dst ~= q{
        auto T()()
        {
            return TransposedAndOrConjugated!(Yes.isTransposed, No.isConjugated, typeof(this))(this);
        }
    };


    if(canFind(list, ".H"))
    dst ~= q{
        static if(isComplex!ElementType)
        {
            auto H()()
            {
                return TransposedAndOrConjugated!(Yes.isTransposed, Yes.isConjugated, typeof(this))(this);
            }
        }
        else
        {
            alias H = T;
        }
    };


    return dst;
}



struct PartialFromMatrix(M, Ridx, Cidx)
if(isMatrixLike!M
    && (is(Ridx == size_t) || is(Ridx == Tuple!(size_t, size_t)))
    && (is(Cidx == size_t) || is(Cidx == Tuple!(size_t, size_t)))
    && (is(Ridx == Tuple!(size_t, size_t)) || is(Cidx == Tuple!(size_t, size_t))))
{
    alias ElementType = M.ElementType;
    enum exprTreeCost = M.exprTreeCost;


    this(M mat, Ridx rs, Cidx cs)
    {
        _mat = mat;
        _rs = rs;
        _cs = cs;

        static if(is(Ridx == Tuple!(size_t, size_t))) {
            _rs[1] = min(_rs[1], mat.length!0);
            assert(_rs[0] < _rs[1]);
        }

        static if(is(Cidx == Tuple!(size_t, size_t))) {
            _cs[1] = min(_cs[1], mat.length!1);
            assert(_cs[0] < _cs[1]);
        }
    }


  static if(is(Ridx == size_t))
  {
    auto ref opIndex(size_t j) inout
    {
        return _mat[_rs, j + _cs[0]];
    }


    size_t length() const { return _cs[1] - _cs[0]; }


    void evalTo(T, SliceKind kind, Alloc)(Slice!(T*, 1, kind) dst, ref Alloc alloc)
    in(dst.length == this.length)
    {
        auto meval = makeViewOrNewSlice(_mat, alloc);
        scope(exit) {
            if(meval.isAllocated)
                alloc.dispose(meval.view.iterator);
        }

        dst[] = meval.view[_rs, _cs[0] .. _cs[1]];
    }


    mixin(definitionsOfVectorOperators(["defaults", "V+V", "M*V", "V*S"]));
  }
  else static if(is(Cidx == size_t))
  {
    auto ref opIndex(size_t i) inout
    {
        return _mat[i + _rs[0], _cs];
    }


    size_t length() const { return _rs[1] - _rs[0]; }


    void evalTo(T, SliceKind kind, Alloc)(Slice!(T*, 1, kind) dst, ref Alloc alloc)
    in(dst.length == this.length)
    {
        auto meval = makeViewOrNewSlice(_mat, alloc);
        scope(exit) {
            if(meval.isAllocated)
                alloc.dispose(meval.view.iterator);
        }

        dst[] = meval.view[_rs[0] .. _rs[1], _cs];
    }


    mixin(definitionsOfVectorOperators(["defaults", "V+V", "M*V", "V*S"]));
  }
  else
  {
    auto ref opIndex(size_t i, size_t j)
    {
        return _mat[i + _rs[0], j + _cs[0]];
    }


    size_t length(size_t dim)() const
    {
        static if(dim == 0)
            return _rs[1] - _rs[0];
        else
            return _cs[1] - _cs[0];
    }


    void evalTo(T, SliceKind kind, Alloc)(Slice!(T*, 2, kind) dst, ref Alloc alloc)
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        auto meval = makeViewOrNewSlice(_mat, alloc);
        scope(exit) {
            if(meval.isAllocated)
                alloc.dispose(meval.view.iterator);
        }

        dst[] = meval.view[_rs[0] .. _rs[1], _cs[0] .. _cs[1]];
    }


    typeof(this) _op_partial_impl_(Tuple!(size_t, size_t) rs1, Tuple!(size_t, size_t) cs1)
    {
        alias TP = Tuple!(size_t, size_t);
        return typeof(this)(_mat, TP(_rs[0] + rs1[0], _rs[0] + rs1[1]), TP(_cs[0] + cs1[0], _cs[0] + cs1[1]));
    }


    auto _op_partial_impl_(Tuple!(size_t, size_t) rs1, size_t j)
    {
        alias TP = Tuple!(size_t, size_t);
        return PartialFromMatrix!(M, TP, size_t)(_mat, TP(_rs[0] + rs1[0], _rs[0] + rs1[1]), _cs[0] + j);
    }


    auto _op_partial_impl_(size_t i, Tuple!(size_t, size_t) cs1)
    {
        alias TP = Tuple!(size_t, size_t);
        return PartialFromMatrix!(M, size_t, TP)(_mat, _rs[0] + i, TP(_cs[0] + cs1[0], _cs[0] + cs1[1]));
    }


    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*V", "M*S", "M*M", ".T", ".H"]));
  }

  private:
    M _mat;
    Ridx _rs;
    Cidx _cs;
}



auto partial(M)(M mat, size_t[2] rs, size_t[2] cs)
if(isMatrixLike!M)
{
    alias TP = Tuple!(size_t, size_t);

    static if(is(typeof(mat._op_partial_impl_(TP(rs[0], rs[1]), TP(cs[0], cs[1])))))
    {
        return mat._op_partial_impl_(TP(rs[0], rs[1]), TP(cs[0], cs[1]));
    }
    else
    {
        return PartialFromMatrix!(M, TP, TP)(mat, TP(rs[0], rs[1]), TP(cs[0], cs[1]));
    }
}


auto partial(M)(M mat, size_t[2] rs, size_t j)
if(isMatrixLike!M)
{
    alias TP = Tuple!(size_t, size_t);

    static if(is(typeof(mat._op_partial_impl_(TP(rs[0], rs[1]), j))))
    {
        return _mat._op_partial_impl_(TP(rs[0], rs[1]), j);
    }
    else
    {
        return PartialFromMatrix!(M, TP, size_t)(_mat, TP(rs[0], rs[1]), j);
    }
}


auto partial(M)(M mat, size_t i, size_t[2] cs)
if(isMatrixLike!M)
{
    alias TP = Tuple!(size_t, size_t);

    static if(is(typeof(mat._op_partial_impl_(i, TP(cs[0], cs[1])))))
    {
        return _mat._op_partial_impl_(i, TP(cs[0], cs[1]));
    }
    else
    {
        return PartialFromMatrix!(M, size_t, TP)(_mat, i, TP(cs[0], cs[1]));
    }
}

unittest
{
    auto mat = matrix!int(3, 3, 0);
    mat[0, 0] = 1;
    mat[0, 1] = 2;
    mat[0, 2] = 3;
    mat[1, 0] = 4;
    mat[1, 1] = 5;
    mat[1, 2] = 6;
    mat[2, 0] = 7;
    mat[2, 1] = 8;
    mat[2, 2] = 9;

    auto p1 = mat.partial([1, 2], [0, 3]);
    assert(p1.length!0 == 1);
    assert(p1.length!1 == 3);
    assert(p1[0, 0] == 4);
    assert(p1[0, 1] == 5);
    assert(p1[0, 2] == 6);

    auto s = slice!int(1, 3);
    p1.evalTo(s, theAllocator);
    assert(s[0, 0] == 4);
    assert(s[0, 1] == 5);
    assert(s[0, 2] == 6);

    auto p2 = p1[0, 1 .. 2];
    assert(p2.length == 1);
    assert(p2[0] == 5);

    auto p3 = p1[0 .. 1, 1];
    assert(p3.length == 1);
    assert(p3[0] == 5);
}
