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


alias matvecAllocator = theAllocator;



/** 
 * 
 */
auto makeViewOrNewSlice(Vec, Alloc)(const Vec vec, ref Alloc alloc)
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


auto makeViewOrNewSlice(Mat, Alloc)(const Mat mat, ref Alloc alloc)
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
}



struct ConstAll(E, E value, size_t Dim)
if(is(typeof(value) : E) && (Dim == 1 || Dim == 2) )
{
    alias ElementType = E;


    this(size_t[Dim] len...){ this._len = len; }


  static if(Dim == 1)
  {
    size_t length() const { return this._len[0]; }

    E opIndex(size_t) const
    {
        return value;
    }

    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
    in(dst.length == this.length)
    {
        dst[] = value;
    }
  }
  else
  {
    size_t length(size_t d)() const { return this._len[d]; }

    E opIndex(size_t, size_t) const
    {
        return value;
    }

    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        dst[] = value;
    }


    auto T() {
        return ConstAll!(E, value, Dim)([_len[1], _len[0]]);
    }
  }


  private:
    size_t[Dim] _len;
}

unittest
{
    auto zeros = ConstAll!(int, 0, 2)();
    static assert(isMatrixLike!(typeof(zeros)));
}


alias ZeroMatrix(E) = ConstAll!(E, E(0), 2);
alias ZeroVector(E) = ConstAll!(E, E(0), 1);


ZeroVector!E zeros(E)(size_t n)
{
    return ZeroVector!E(n);
}


ZeroMatrix!E zeros(E)(size_t n, size_t m)
{
    return ZeroMatrix!E(n, m);
}


enum isZeros(T) = is(Unqual!T == ZeroMatrix!(T.ElementType)) || is(Unqual!T == ZeroVector!(T.ElementType));


struct ConstEye(E, E value)
{
    alias ElementType = E;


    this(size_t len){ this._len = len; }


    size_t length(size_t d)() const { return _len; }


    E opIndex(size_t i, size_t j) const
    {
        if(i == j)
            return value;
        else
            return E(0);
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        dst[] = value;
    }


    auto T() {
        return this;
    }


  private:
    size_t _len;
}


unittest
{
    auto ident = ConstEye!(int, 3)(2);
    static assert(isMatrixLike!(typeof(ident)));

    assert(ident[0, 0] == 3);
    assert(ident[0, 1] == 0);
    assert(ident[1, 0] == 0);
    assert(ident[1, 1] == 3);
}


alias Identity(E) = ConstEye!(E, E(1));


Identity!E identity(E)(size_t n)
{
    return Identity!E(n);
}


enum isIdentity(T) = is(Unqual!T == Identity!(T.ElementType));



// 
struct MatrixMatrixMulGEMM(S, MatA, MatB, MatC)
{
    alias ElementType = typeof(S.init * MatA.ElementType.init * MatB.ElementType.init + S.init * MatC.ElementType.init);


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


    mixin MatrixOperators!(["M*V", "M+M"]) OpImpls;
    mixin(definitionsOfMatrixOperators(["M*M", ".H", ".T"]));


    auto opBinary(string op : "*", U)(U u)
    if(!isMatrixLike!U && !isVectorLike!U && is(typeof(U.init * S.init)))
    {
        return matrixGemm(_alpha * u, _matA, _matB, _beta * u, _matC);
    }


    auto opBinaryRight(string op : "*", U)(U u)
    if(!isMatrixLike!U && !isVectorLike!U && is(typeof(U.init * S.init)))
    {
        return matrixGemm(_alpha * u, _matA, _matB, _beta * u, _matC);
    }


    auto opBinary(string op : "*", V)(V vec)
    if(isVectorLike!V)
    in(vec.length == this.length!1)
    {
        return _matA * (_matB * vec) * _alpha + (_matC * vec) * _beta;
    }


    auto opBinary(string op : "+", M)(M mat)
    if(isMatrixLike!M)
    in(mat.length!0 == this.length!0 && mat.length!1 == this.length!1)
    {
        static if(isZeros!MatC)
        {
            static if(is(M == MatrixAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
                return matrixGemm(_alpha, _matA, _matB, mat._alpha, mat._matA);
            else
                return matrixGemm(_alpha, _matA, _matB, ElementType(1), mat);
        }
        else
            return OpImpls.opBinary!op(mat);
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
if(isMatrixLike!MatA && isVectorLike!VecB)
{
    alias ElementType = typeof(T.init * MatA.ElementType.init * VecB.ElementType.init + T.init * VecC.ElementType.init);


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


    auto opBinary(string op : "*", S)(S val)
    if(is(T :  ElementType))
    {
        return vectorGemv(_alpha * val, _matA, _vecB, _beta * val, _vecC);
    }


    auto opBinaryRight(string op : "*", S)(S val)
    if(is(T : ElementType))
    {
        return vectorGemv(_alpha * val, _matA, _vecB, _beta * val, _vecC);
    }


    static if(isZeros!VecC)
    {
        auto opBinary(string op : "+", V)(V vec)
        if(isVectorLike!V)
        in(vec.length == this.length)
        {
            // pragma(msg, V);
            static if(is(V == VectorAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
                return vectorGemv(_alpha, _matA, _vecB, vec._alpha, vec._vecA);
            else
                return vectorGemv(_alpha, _matA, _vecB, ElementType(1), vec);
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
        auto viewA = _matA.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        _matB.evalTo(dst, alloc);
        dst[] *= _beta;

        foreach(i; 0 .. this.length!0)
            foreach(j; 0 .. this.length!1)
                dst[i, j] += _alpha * viewA.view[i, j];
    }


    mixin MatrixOperators!(["M+M",]) OpImpls;
    mixin(definitionsOfMatrixOperators(["M*M", "M*V", ".H", ".T"]));


    auto opBinary(string op : "*", U)(U u)
    if(!isVectorLike!U && !isMatrixLike!U && is(typeof(U.init * S.init)))
    {
        return matrixAxpby(_alpha * u, _matA, _beta * u, _matB);
    }


    auto opBinaryRight(string op : "*", U)(U u)
    if(!isVectorLike!U && !isMatrixLike!U && is(typeof(U.init * S.init)))
    {
        return matrixAxpby(_alpha * u, _matA, _beta * u, _matB);
    }


    auto opBinary(string op : "+", M)(M mat)
    if(isMatrixLike!M)
    in(mat.length!0 == this.length!0 && mat.length!1 == this.length!1)
    {
        static if(isZeros!MatB)
        {
            return matrixAxpby(_alpha, _matA, S(0), mat);
        }
        else
            return OpImpls.opBinary!op(mat);
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
        auto viewA = _vecA.makeViewOrNewSlice(alloc);
        scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

        _vecB.evalTo(dst, alloc);
        dst[] *= _beta;

        foreach(i; 0 .. this.length)
            dst[i] += _alpha * viewA.view[i];
    }


    mixin VectorOperators!(["V+V"]) OpImpls;


    auto opBinary(string op : "*", U)(U u)
    if(!isVectorLike!U && !isMatrixLike!U && is(typeof(U.init * T.init)))
    {
        return vectorAxpby(_alpha * u, _matA, _beta * u, _matB);
    }


    auto opBinaryRight(string op : "*", U)(U u)
    if(!isVectorLike!U && !isMatrixLike!U && is(typeof(U.init * T.init)))
    {
        return vectorAxpby(_alpha * u, _matA, _beta * u, _matB);
    }


    static if(isZeros!VecB)
    {
        auto opBinary(string op : "+", V)(V vec)
        if(isVectorLike!V)
        in(vec.length == vec.length)
        {
            static if(is(V == VectorAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
                return vectorAxpby(_alpha, _vecA, vec._alpha, vec._vecA);
            else
                return vectorAxpby(_alpha, _vecA, T(1), vec);
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


    mixin MatrixOperators!(["M*V", "M+M", ".T", ".H"]);
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


    mixin VectorOperators!(["S*V", "V*S", "V+V"]);
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


    mixin MatrixOperators!(["M*V", "M+M", ".H"]);


  private:
    Mat _mat;
}


Transposed!Mat matrixTransposed(Mat)(Mat mat)
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
    };

    if(canFind(list, "S*V"))
    dst ~= q{
        auto opBinary(string op : "*", S)(S scalar)
        if(!isVectorLike!S && !isMatrixLike!S && is(typeof(S.init * ElementType.init)))
        {
            return vectorAxpby(scalar, this, S(0), zeros!S(this.length));
        }
    };


    if(canFind(list, "V*S"))
    dst ~= q{
        auto opBinaryRight(string op : "*", S)(S scalar)
        if(!isVectorLike!S && !isMatrixLike!S && is(typeof(S.init * ElementType.init)))
        {
            return vectorAxpby(scalar, this, S(0), zeros!T(this.length));
        }
    };


    if(canFind(list, "V+V"))
    dst ~= q{
        auto opBinary(string op : "+", V)(V vec)
        if(isVectorLike!V)
        in(vec.length == this.length)
        {
            static if(isZeros!V)
                return this;
            else
                return vectorAxpby(ElementType(1), this, ElementType(1), vec);
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
    };


    if(canFind(list, "M*M"))
    dst ~= q{
        auto opBinary(string op : "*", M)(M mat)
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


    if(canFind(list, "M*S"))
    dst ~= q{
        auto opBinary(string op : "*", S)(S scalar)
        if(!isVectorLike!S && !isMatrixLike!S && is(typeof(S.init * ElementType.init)))
        {
            return matrixAxpby(scalar, this, S(0), zeros!ElementType(this.length!0, this.length!1));
        }
    };


    if(canFind(list, "M+M"))
    dst ~= q{
        auto opBinary(string op : "+", M)(M mat)
        if(isMatrixLike!M)
        in(mat.length!0 == this.length!0 && mat.length!1 == this.length!1)
        {
            static if(isZeros!M)
                return this;
            else {
                static if(is(M == MatrixAxpby!(A, V1, V2), A, V1, V2) && isZeros!V2)
                    return matrixAxpby(ElementType(1), this, mat._alpha, mat._matA);
                else
                    return matrixAxpby(ElementType(1), this, S(1), mat);
            }
        }
    };


    if(canFind(list, "M*V"))
    dst ~= q{
        auto opBinary(string op : "*", V)(V vec)
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
