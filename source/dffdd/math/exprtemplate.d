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

enum hasMemoryView(T) = (is(typeof(T.init.sliced().lightScope())) && hasMemoryView!(typeof(T.init.sliced().lightScope())))
                     || is(Unqual!(typeof(T.init.sliced())) == Slice!(T.ElementType*, 1, kind), SliceKind kind)
                     || (is(Unqual!(typeof(T.init.sliced())) == Slice!(T.ElementType*, 2, kind), SliceKind kind) && (kind == Contiguous || kind == Canonical));

/** 
 * 
 */
auto makeViewOrNewSlice(Vec, Alloc)(const Vec vec, ref Alloc alloc)
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


auto makeViewOrNewSlice(Mat, Alloc)(const Mat mat, ref Alloc alloc)
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


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
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


    mixin(definitionsOfMatrixOperators(["defaults", "M*M", ".H", ".T"]));


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
        return _matA * (_matB * vec) * _alpha + (_matC * vec) * _beta;
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
  else
  {
    mixin(definitionsOfMatrixOperators(["M+M"]));
  }


  private:
    S _alpha;
    S _beta;
    MatA _matA;
    MatB _matB;
    MatC _matC;
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

    auto z1 = mat1 * mat2 * 2 + mat3 * 3;
    assert(z1.makeViewOrNewSlice(theAllocator).view == [[20, 29], [62, 89]]);

    auto z2 = z1 * vec1;
    assert(z2.makeViewOrNewSlice(theAllocator).view == [127, 391]);
}




// 
struct MatrixVectorMulGEMV(T, MatA, VecB, VecC)
if(isMatrixLike!MatA && isVectorLike!VecB && !isMatrixLike!T && !isVectorLike!T)
{
    alias ElementType = typeof(T.init * MatA.ElementType.init * VecB.ElementType.init + T.init * VecC.ElementType.init);
    // pragma(msg, "ElementType is " ~ ElementType.stringof);
    // pragma(msg, "T is " ~ T.stringof);
    // pragma(msg, "MatA.ElementType is " ~ MatA.ElementType.stringof);
    // pragma(msg, "VecB.ElementType is " ~ VecB.ElementType.stringof);
    // pragma(msg, "VecC.ElementType is " ~ VecC.ElementType.stringof);
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


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
    in(dst.length == this.length)
    {
        auto viewA = _matA.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        auto viewB = _vecB.makeViewOrNewSlice(alloc);
        scope(exit) if(viewB.isAllocated) alloc.dispose(viewB.view.iterator);

        _vecC.evalTo(dst, alloc);

        static if(isBlasType!T)
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
        if(isVectorLike!V)
        in(vec.length == this.length)
        {
            static if(is(V == VectorAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
                return vectorGemv(_alpha, _matA, _vecB, op == "+" ? vec._alpha : -vec._alpha, vec._vecA);
            else
                return vectorGemv(_alpha, _matA, _vecB, ElementType(op == "+" ? 1 : -1), vec);
        }

        auto _opBR_Impl_(string op, V)(V vec)
        if(isVectorLike!V)
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


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
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
    mixin(definitionsOfMatrixOperators(["defaults", ".H", ".T"]));


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
  else
  {
    mixin(definitionsOfMatrixOperators(["M+M", "M*M", "M*V"]));
  }


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


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
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

            foreach(i; 0 .. this.length)
                dst[i] += _alpha * viewA.view[i];
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


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
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


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
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


// transpose
struct Transposed(Mat)
if(isMatrixLike!Mat)
{
    alias ElementType = Mat.ElementType;
    enum float exprTreeCost = Mat.exprTreeCost + 1;

    this(Mat mat)
    {
        _mat = mat;
    }


    ElementType opIndex(size_t i, size_t j) const
    {
        return _mat[j, i];
    }


    size_t length(size_t dim)() const
    if(dim == 0 || dim == 1)
    {
        static if(dim == 0)
            return _mat.length!1;
        else
            return _mat.length!0;
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        auto viewA = _mat.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        dst[] = viewA.view.transposed;
    }


    Mat T()
    {
        return _mat;
    }


    mixin MatrixOperators!(["defaults", "M*M", "M*V", "M+M", ".H"]);


  private:
    Mat _mat;
}


Transposed!Mat matrixTransposed(Mat)(in Mat mat)
{
    return Transposed!Mat(mat);
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
        if(!isMatrixLike!X && !isVectorLike!X && is(X : ElementType))
        {
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
    import std.traits : isFloatingPoint, isIntegral;
    import dffdd.math.linalg : isBlasType;
    import dffdd.math.exprtemplate;
    import dffdd.math.matrix;
    import dffdd.math.vector;
    import dffdd.math.matrixspecial;
    };


    if(canFind(list, "defaults"))
    dst ~= q{
        auto opBinary(string op, X)(X rhs)
        if(op == "+" || op == "-" || op == "*")
        {
            static if(is(typeof(this._opB_Impl_!op(rhs))) && is(typeof(rhs._opBR_Impl_!op(this))))
            {
                pragma(msg, typeof(this));
                pragma(msg, X);
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
        if(!isMatrixLike!X && !isVectorLike!X && is(X : ElementType))
        {
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
                    return matrixAxpby(ElementType(1), this, S(op == "+" ? 1 : -1), mat);
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
            return matrixTransposed(this);
        }
    };


    if(canFind(list, ".H"))
    dst ~= q{
        auto H()()
        {
            import dffdd.math.complex : isComplex, conj;
            static if(isComplex!ElementType) {
                return elementMap!conj(this.T());
            } else {
                return this.T();
            }
        }
    };


    return dst;
}
