module dffdd.math.matrixspecial;

import std.algorithm : max, min;
import std.experimental.allocator;
import std.typecons;
import std.traits;
import std.meta;

import mir.ndslice;

import dffdd.math.vector;
import dffdd.math.matrix;
import dffdd.math.exprtemplate;
import dffdd.math.complex;
import dffdd.math.linalg : isBlasType;


/+
auto allSameVector(E)(size_t n, E value)
{
    return AllSameElements!(E, 1)([n], value);
}


auto allSameMatrix(E)(size_t n, size_t m, E value)
{
    return AllSameElements!(E, 2)([n, m], value);
}


struct AllSameElements(E, size_t Dim)
if(is(typeof(value) : E) && (Dim == 1 || Dim == 2))
{
    alias ElementType = E;
    enum float exprTreeCost = 1;

    this(size_t[Dim] len, E value){
        this._len = len;
        this._value = value;
    }


  static if(Dim == 1)
  {
    size_t length() const { return this._len[0]; }

    E opIndex(size_t) const
    {
        return _value;
    }

    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc)
    in(dst.length == this.length)
    {
        dst[] = _value;
    }
  }
  else
  {
    size_t length(size_t d)() const { return this._len[d]; }

    E opIndex(size_t, size_t) const
    {
        return _value;
    }

    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc)
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        dst[] = _value;
    }

    auto T() {
        return AllSameElements!(E, Dim)([_len[1], _len[0]], _value);
    }

    auto H() {
      static if(is(typeof(E.init.re)))
      {
        return AllSameElements!(E, Dim)([_len[1], _len[0]], conj(_value));
      }
      else
      {
         return AllSameElements!(E, Dim)([_len[1], _len[0]], value); 
      }
    }
  }

  static if(Dim == 1)
  {
    mixin(definitionsOfVectorOperators(["defaults"]));
    mixin VectorOperators!(["V+V"]) OpImpls;

    auto _opB_Impl_(string op, V)(V vec)
    if(isVectorLike!V && (op == "+" || op == "-"))
    in(this._length == vec._length)
    {
        static if(is(V == AllSameElements!(E1, 1), E1)) {
            return allSameVector(this._length, this._value + vec._value);
        } else {
            return OpImpls._opB_Impl_(vec);
        }
    }

    auto _opB_Impl_(string op : "*", S)(S scalar)
    if(!isMatrixLike!S && !isVectorLike!S)
    {
        return allSameVector(this._length, this._value * scalar);
    }
  }
  else
  {
    mixin(definitionsOfMatrixOperators(["defaults", "M*V"]));
    mixin MatrixOperators!(["M+M"]) OpImpls;

    auto _opB_Impl_(string op, V)(V vec)
    if(isVectorLike!V && (op == "+" || op == "-"))
    in(this._length == vec._length)
    {
        static if(is(V == AllSameElements!(E1, 2), E1)) {
            return allSameMatrix(this._length, this._value + vec._value);
        } else {
            return OpImpls._opB_Impl_(vec);
        }
    }

    auto _opB_Impl_(string op : "*", S)(S scalar)
    if(!isMatrixLike!S && !isVectorLike!S)
    {
        return allSameMatrix(this._length, this._value * scalar);
    }
  }


  private:
    size_t[Dim] _length;
    E _value;
}
+/


struct ConstAll(E, E value, size_t Dim)
if(is(typeof(value) : E) && (Dim == 1 || Dim == 2) )
{
    alias ElementType = E;
    enum float exprTreeCost = 1;


    this(size_t[Dim] len...){ this._len = len; }


  static if(Dim == 1)
  {
    size_t length() const { return this._len[0]; }

    E opIndex(size_t) const
    {
        return value;
    }

    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc)
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

    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc)
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        dst[] = value;
    }


    auto T() {
        return ConstAll!(E, value, Dim)([_len[1], _len[0]]);
    }
  }

  static if(Dim == 1)
  {
    mixin(definitionsOfVectorOperators(["defaults", "V+V", "V*S"]));
  }
  else
  {
    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*M", "M*V", "M*S"]));
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


ConstAll!(E, E(1), 1) ones(E)(size_t n)
{
    return ConstAll!(E, E(1), 1)(n);
}


ConstAll!(E, E(1), 2) ones(E)(size_t n, size_t m)
{
    return ConstAll!(E, E(1), 2)(n, m);
}


enum isZeros(T) = is(Unqual!T == ZeroMatrix!(T.ElementType)) || is(Unqual!T == ZeroVector!(T.ElementType));


struct ConstEye(E, E value)
{
    alias ElementType = E;
    enum float exprTreeCost = 1;


    this(size_t len){ this._len = len; }


    size_t length(size_t d)() const { return _len; }


    E opIndex(size_t i, size_t j) const
    {
        if(i == j)
            return value;
        else
            return E(0);
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc)
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        import mir.ndslice : diagonal;

        dst[] = T(0);
        dst.diagonal[] = cast(T)value;
    }


    auto T() {
        return this;
    }


    static if(value == 1)
    {
        auto _opB_Impl_(string op : "*", V)(V vec)
        if(isVectorLike!V)
        in(vec.length == this.length!1)
        {
            return vec;
        }

        auto _opB_Impl_(string op : "*", M)(M mat)
        if(isMatrixLike!M)
        in(mat.length!1 == this.length!0)
        {
            return mat;
        }

        auto _opBR_Impl_(string op : "*", M)(M mat)
        if(isMatrixLike!M)
        in(mat.length!0 == this.length!1)
        {
            return mat;
        }
    }
    else
    {
        auto _opB_Impl_(string op : "*", V)(V vec)
        if(isVectorLike!V)
        in(vec.length == this.length!1)
        {
            return vec * value;
        }

        auto _opB_Impl_(string op : "*", M)(M mat)
        if(isMatrixLike!M)
        in(mat.length!1 == this.length!0)
        {
            return mat * value;
        }

        auto _opBR_Impl_(string op : "*", M)(M mat)
        if(isMatrixLike!M)
        in(mat.length!0 == this.length!1)
        {
            return mat * value;
        }
    }


    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*S"]));


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

unittest
{
    import std.complex;
    auto ident = ConstEye!(Complex!float, Complex!float(3, 0))(2);
    static assert(isMatrixLike!(typeof(ident)));

    assert(ident[0, 0].re == 3);
    assert(ident[0, 1].re == 0);
    assert(ident[1, 0].re == 0);
    assert(ident[1, 1].re == 3);
}


alias Identity(E) = ConstEye!(E, E(1));


Identity!(Unqual!E) identity(E)(size_t n)
{
    return Identity!(Unqual!E)(n);
}


enum isIdentity(T) = is(Unqual!T == Identity!(T.ElementType));


struct MulDiagonal(ArrayLike, MatVec)
if(is(typeof((ArrayLike a){ auto e1 = a[0]; auto e2 = a.length; })) && (isMatrixLike!MatVec || isVectorLike!MatVec))
{
    alias ElementType = typeof(_diag[0] * MatVec.ElementType.init);
    enum float exprTreeCost = MatVec.exprTreeCost + EXPR_COST_N;


    this(size_t M, size_t N, ArrayLike vec, MatVec matvec)
    // in((vec.length == M || vec.length == N) && matvec.length!0 == N)
    in {
        assert(vec.length == min(M, N));
        static if(isMatrixLike!MatVec)
            assert(matvec.length!0 == N);
        else
            assert(matvec.length == N);
    }
    do {
        _M = M;
        _N = N;
        _diag = vec;
        _matvec = matvec;
    }


  static if(isMatrixLike!MatVec)
  {
    size_t length(size_t dim = 0)() const @property
    if(dim == 0 || dim == 1)
    {
        static if(dim == 0)
            return _M;
        else
            return _matvec.length!1;
    }


    ElementType opIndex(size_t i, size_t j) const
    {
        if(i < _diag.length)
            return _matvec[i, j] * _diag[i];
        else
            return ElementType(0);
    }
  }
  else
  {
    size_t length() const @property
    {
        // return _matvec.length;
        return _M;
    }


    ElementType opIndex(size_t i) const
    {
        if(i < _diag.length)
            return _matvec[i] * _diag[i];
        else
            return ElementType(0);
    }
  }


  static if(isMatrixLike!MatVec)
  {
    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc)
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        if(_M == _N) {
            _matvec.evalTo(dst, alloc);
            foreach(i; 0 .. _diag.length)
                dst[i, 0 .. $] *= _diag[i];
        } else {
            auto viewA = _matvec.makeViewOrNewSlice(alloc);
            scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

            dst[] = T(0);
            foreach(i; 0 .. _diag.length)
                dst[i, 0 .. $] = _diag[i] * viewA.view[i, 0 .. $];
        }
    }
  }
  else
  {
    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc)
    in(dst.length == this.length)
    {
        if(_M == _N) {
            _matvec.evalTo(dst, alloc);
            foreach(i; 0 .. _diag.length)
                dst[i] *= _diag[i];
        } else {
            auto viewA = _matvec.makeViewOrNewSlice(alloc);
            scope(exit) if(viewA.isAllocated) alloc.dispose(viewA.view.iterator);

            dst[] = T(0);
            foreach(i; 0 .. _diag.length)
                dst[i] = _diag[i] * viewA.view[i];
        }
    }
  }


  static if(isMatrixLike!MatVec)
  {
    auto _opB_Impl_(string op : "*", V)(V vec)
    if(isVectorLike!V)
    in(vec.length == this.length!1)
    {
        auto mm = _matvec * vec;
        return MulDiagonal!(ArrayLike, typeof(mm))(_M, _N, _diag, mm);
    }


    auto _opB_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(mat.length!0 == this.length!1)
    {
        auto mm = _matvec * mat;
        return MulDiagonal!(ArrayLike, typeof(mm))(_M, _N, _diag, mm);
    }

    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*S"]));
  }
  else
  {
    mixin(definitionsOfVectorOperators(["defaults", "V+V", "V*S"]));
  }


  private:
    size_t _M, _N;
    ArrayLike _diag;
    MatVec _matvec;
}


auto diag(Iterator, SliceKind kind)(in Slice!(Iterator, 1, kind) slice)
{
    return .diag(slice.length, slice.length, slice);
}


auto diag(Iterator, SliceKind kind)(size_t M, size_t N, in Slice!(Iterator, 1, kind) slice)
in(slice.length == min(M, N))
{
    auto s = slice.lightConst;
    auto ident = identity!(typeof(slice[0]))(N);
    return MulDiagonal!(typeof(s), typeof(ident))(M, N, s, ident);
}


auto diag(E)(in E[] slice)
{
    return .diag(slice.length, slice.length, slice);
}


auto diag(E)(size_t M, size_t N, in E[] slice)
in(slice.length == min(M, N))
{
    auto ident = identity!(Unqual!(typeof(slice[0])))(N);
    return MulDiagonal!(const(E)[], typeof(ident))(M, N, slice, ident);
}


unittest
{
    auto ident = identity!int(3);
    static assert(isMatrixLike!(typeof(ident)));

    auto ms = slice!int(3, 3);

    auto s1 = sliced([1, 2, 3]);
    auto dmat1 = diag(s1);

    assert(dmat1.length!0 == 3 && dmat1.length!1 == 3);
    assert(dmat1[0, 0] == 1);
    assert(dmat1[0, 1] == 0);
    assert(dmat1[0, 2] == 0);
    assert(dmat1[1, 0] == 0);
    assert(dmat1[1, 1] == 2);
    assert(dmat1[1, 2] == 0);
    assert(dmat1[2, 0] == 0);
    assert(dmat1[2, 1] == 0);
    assert(dmat1[2, 2] == 3);

    dmat1.evalTo(ms, theAllocator);
    assert(ms[0, 0] == 1);
    assert(ms[0, 1] == 0);
    assert(ms[0, 2] == 0);
    assert(ms[1, 0] == 0);
    assert(ms[1, 1] == 2);
    assert(ms[1, 2] == 0);
    assert(ms[2, 0] == 0);
    assert(ms[2, 1] == 0);
    assert(ms[2, 2] == 3);

    auto s2 = [1, 2, 3];
    auto dmat2 = diag(s2);

    auto mat3 = [1, 2, 3, 4, 5, 6, 7, 8, 9].sliced(3, 3).matrixed;
    auto mul3 = dmat2 * mat3;
    assert(mul3[0, 0] == 1);
    assert(mul3[0, 1] == 2);
    assert(mul3[0, 2] == 3);
    assert(mul3[1, 0] == 8);
    assert(mul3[1, 1] == 10);
    assert(mul3[1, 2] == 12);
    assert(mul3[2, 0] == 21);
    assert(mul3[2, 1] == 24);
    assert(mul3[2, 2] == 27);

    mul3.evalTo(ms, theAllocator);
    assert(ms[0, 0] == 1);
    assert(ms[0, 1] == 2);
    assert(ms[0, 2] == 3);
    assert(ms[1, 0] == 8);
    assert(ms[1, 1] == 10);
    assert(ms[1, 2] == 12);
    assert(ms[2, 0] == 21);
    assert(ms[2, 1] == 24);
    assert(ms[2, 2] == 27);

    auto v1 = vector!int(3);
    // auto mul4 = mul3 * v1;
    auto mul4 = mul3._opB_Impl_!"*"(v1);
}

unittest
{
    auto mat1 = [1, 2, 3, 4, 5, 6, 7, 8, 9].sliced(3, 3).matrixed;
    auto diag1 = diag(2, 3, [1, 2]);
    auto mul1 = diag1 * mat1;
    assert(mul1.length!0 == 2);
    assert(mul1.length!1 == 3);
    assert(mul1[0, 0] == 1);
    assert(mul1[0, 1] == 2);
    assert(mul1[0, 2] == 3);
    assert(mul1[1, 0] == 8);
    assert(mul1[1, 1] == 10);
    assert(mul1[1, 2] == 12);
    
    auto s1 = slice!int(2, 3);
    mul1.evalTo(s1, theAllocator);
    assert(s1[0, 0] == 1);
    assert(s1[0, 1] == 2);
    assert(s1[0, 2] == 3);
    assert(s1[1, 0] == 8);
    assert(s1[1, 1] == 10);
    assert(s1[1, 2] == 12);

    auto diag2 = diag(4, 3, [2, 1, 2]);
    auto mul2 = diag2 * mat1;
    assert(mul2.length!0 == 4);
    assert(mul2.length!1 == 3);
    assert(mul2[0, 0] == 2);
    assert(mul2[0, 1] == 4);
    assert(mul2[0, 2] == 6);
    assert(mul2[1, 0] == 4);
    assert(mul2[1, 1] == 5);
    assert(mul2[1, 2] == 6);
    assert(mul2[2, 0] == 14);
    assert(mul2[2, 1] == 16);
    assert(mul2[2, 2] == 18);
    assert(mul2[3, 0] == 0);
    assert(mul2[3, 1] == 0);
    assert(mul2[3, 2] == 0);

    auto s2 = slice!int(4, 3);
    mul2.evalTo(s2, theAllocator);
    assert(s2[0, 0] == 2);
    assert(s2[0, 1] == 4);
    assert(s2[0, 2] == 6);
    assert(s2[1, 0] == 4);
    assert(s2[1, 1] == 5);
    assert(s2[1, 2] == 6);
    assert(s2[2, 0] == 14);
    assert(s2[2, 1] == 16);
    assert(s2[2, 2] == 18);
    assert(s2[3, 0] == 0);
    assert(s2[3, 1] == 0);
    assert(s2[3, 2] == 0);
}


/**
Unitary DFTMatrix
*/
struct DFTMatrix(C, Flag!"isFwd" isFwd = Yes.isFwd)
if(isNarrowComplex!C)
{
    alias ElementType = C;
    enum float exprTreeCost = 0;

    import std.traits : TemplateOf;
    private alias CpxTemplate = TemplateOf!C;

    import dffdd.utils.fft;

    /**
    make Unitary DFT Matrix
    */
    this(size_t n)
    {
        _N = n;
        _fftw = makeFFTWObject!CpxTemplate(_N);
    }


    size_t length(size_t dim)() const if(dim == 0 || dim == 1) { return _N; }


    auto opIndex(size_t i, size_t j) const
    {
        import std.math : PI;
        import dffdd.math.math : fast_expi, fast_sqrt;

        return fast_expi!C(2*PI/_N*i*j * (isFwd ? -1 : 1)) / fast_sqrt!(typeof(C.init.re))(_N*1.0);
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        foreach(i; 0 .. _N)
            foreach(j; 0 .. _N)
                dst[i, j] = this[i, j];
    }


    auto H()
    {
        DFTMatrix!(C, isFwd ? No.isFwd : Yes.isFwd) dst;
        dst._N = this._N;
        dst._fftw = this._fftw;
        return dst;
    }


    auto H() const
    {
        return DFTMatrix!(C, isFwd ? No.isFwd : Yes.isFwd)(_N);
    }


    auto lightConst() const
    {
        return DFTMatrix!(C, isFwd)(_N);
    }


    // auto _opB_Impl_(string op : "*", V)(V vec)
    // if(isVectorLike!V)
    // in(vec.length == this.length!1)
    // {
    //     auto mm = _matvec * vec;
    //     return MulDiagonal!(ArrayLike, typeof(mm))(_M, _N, _diag, mm);
    // }


    // auto _opB_Impl_(string op : "*", M)(M mat)
    // if(isMatrixLike!M)
    // in(mat.length!0 == this.length!1)
    // {
    //     auto mm = _matvec * mat;
    //     return MulDiagonal!(ArrayLike, typeof(mm))(_M, _N, _diag, mm);
    // }


    import dffdd.math.exprtemplate;
    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*V", "M*S", ".T"]));


  private:
    size_t _N;
    FFTWObject!CpxTemplate _fftw;
}


auto dftMatrix(C)(size_t n)
{
    return DFTMatrix!(C, Yes.isFwd)(n);
}


auto idftMatrix(C)(size_t n)
{
    return DFTMatrix!(C, No.isFwd)(n);
}


unittest
{
    import std.math : isClose;

    alias C = Complex!double;
    auto dftmat = dftMatrix!C(4) * 2;
    auto idftmat = dftmat.H;

    enum float REPS = 0,        // 相対誤差でのチェックを無効
               AEPS = 1e-9;     // 絶対誤差でチェック

    assert(isClose(dftmat[0, 0].re, 1, REPS, AEPS) && isClose(dftmat[0, 0].im, 0, REPS, AEPS));
    assert(isClose(dftmat[0, 1].re, 1, REPS, AEPS) && isClose(dftmat[0, 1].im, 0, REPS, AEPS));
    assert(isClose(dftmat[0, 2].re, 1, REPS, AEPS) && isClose(dftmat[0, 2].im, 0, REPS, AEPS));
    assert(isClose(dftmat[0, 3].re, 1, REPS, AEPS) && isClose(dftmat[0, 3].im, 0, REPS, AEPS));

    assert(isClose(dftmat[1, 0].re, 1, REPS, AEPS) && isClose(dftmat[1, 0].im, 0, REPS, AEPS));
    assert(isClose(dftmat[1, 1].re, 0, REPS, AEPS) && isClose(dftmat[1, 1].im, -1, REPS, AEPS));
    assert(isClose(dftmat[1, 2].re, -1, REPS, AEPS) && isClose(dftmat[1, 2].im, 0, REPS, AEPS));
    assert(isClose(dftmat[1, 3].re, 0, REPS, AEPS) && isClose(dftmat[1, 3].im, 1, REPS, AEPS));

    assert(isClose(dftmat[2, 0].re, 1, REPS, AEPS) && isClose(dftmat[2, 0].im, 0, REPS, AEPS));
    assert(isClose(dftmat[2, 1].re, -1, REPS, AEPS) && isClose(dftmat[2, 1].im, 0, REPS, AEPS));
    assert(isClose(dftmat[2, 2].re, 1, REPS, AEPS) && isClose(dftmat[2, 2].im, 0, REPS, AEPS));
    assert(isClose(dftmat[2, 3].re, -1, REPS, AEPS) && isClose(dftmat[2, 3].im, 0, REPS, AEPS));

    assert(isClose(dftmat[3, 0].re, 1, REPS, AEPS) && isClose(dftmat[3, 0].im, 0, REPS, AEPS));
    assert(isClose(dftmat[3, 1].re, 0, REPS, AEPS) && isClose(dftmat[3, 1].im, 1, REPS, AEPS));
    assert(isClose(dftmat[3, 2].re, -1, REPS, AEPS) && isClose(dftmat[3, 2].im, 0, REPS, AEPS));
    assert(isClose(dftmat[3, 3].re, 0, REPS, AEPS) && isClose(dftmat[3, 3].im, -1, REPS, AEPS));

    assert(isClose(idftmat[0, 0].re, 1, REPS, AEPS) && isClose(idftmat[0, 0].im, 0, REPS, AEPS));
    assert(isClose(idftmat[0, 1].re, 1, REPS, AEPS) && isClose(idftmat[0, 1].im, 0, REPS, AEPS));
    assert(isClose(idftmat[0, 2].re, 1, REPS, AEPS) && isClose(idftmat[0, 2].im, 0, REPS, AEPS));
    assert(isClose(idftmat[0, 3].re, 1, REPS, AEPS) && isClose(idftmat[0, 3].im, 0, REPS, AEPS));

    assert(isClose(idftmat[1, 0].re, 1, REPS, AEPS) && isClose(idftmat[1, 0].im, 0, REPS, AEPS));
    assert(isClose(idftmat[1, 1].re, 0, REPS, AEPS) && isClose(idftmat[1, 1].im, 1, REPS, AEPS));
    assert(isClose(idftmat[1, 2].re, -1, REPS, AEPS) && isClose(idftmat[1, 2].im, 0, REPS, AEPS));
    assert(isClose(idftmat[1, 3].re, 0, REPS, AEPS) && isClose(idftmat[1, 3].im, -1, REPS, AEPS));

    assert(isClose(idftmat[2, 0].re, 1, REPS, AEPS) && isClose(idftmat[2, 0].im, 0, REPS, AEPS));
    assert(isClose(idftmat[2, 1].re, -1, REPS, AEPS) && isClose(idftmat[2, 1].im, 0, REPS, AEPS));
    assert(isClose(idftmat[2, 2].re, 1, REPS, AEPS) && isClose(idftmat[2, 2].im, 0, REPS, AEPS));
    assert(isClose(idftmat[2, 3].re, -1, REPS, AEPS) && isClose(idftmat[2, 3].im, 0, REPS, AEPS));

    assert(isClose(idftmat[3, 0].re, 1, REPS, AEPS) && isClose(idftmat[3, 0].im, 0, REPS, AEPS));
    assert(isClose(idftmat[3, 1].re, 0, REPS, AEPS) && isClose(idftmat[3, 1].im, -1, REPS, AEPS));
    assert(isClose(idftmat[3, 2].re, -1, REPS, AEPS) && isClose(idftmat[3, 2].im, 0, REPS, AEPS));
    assert(isClose(idftmat[3, 3].re, 0, REPS, AEPS) && isClose(idftmat[3, 3].im, 1, REPS, AEPS));
}


// struct MulDFTMatrix