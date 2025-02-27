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


  static if(isComplex!E)
  {
    auto H() {
        return ConstEye!(E, conj(value))(_len);
    }
  }
  else
  {
    alias H = T;
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

    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*S", ".T", ".H"]));
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
    enum float exprTreeCost = 1;

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


    auto T() inout
    {
        return this;
    }


    auto lightConst() const
    {
        return DFTMatrix!(C, isFwd)(_N);
    }


    auto _opB_Impl_(string op : "*", V)(V vec)
    if(isVectorLike!V)
    in(vec.length == this.length!1)
    {
        return DFTMatrixMulMV!(C, isFwd, V)(this, vec);
    }


    auto _opB_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(this.length!1 == mat.length!0)
    {
        static if(is(M == .DFTMatrix!(C, isFwd ? No.isFwd : Yes.isFwd)))
            return identity!C(_N);
        else
            return DFTMatrixMulMM!(C, isFwd, Yes.applyFromLHS, M)(this, mat);
    }


    auto _opBR_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(this.length!0 == mat.length!1)
    {
        return DFTMatrixMulMM!(C, isFwd, No.applyFromLHS, M)(this, mat);
    }


    import dffdd.math.exprtemplate;
    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*S"]));


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

    static assert(isMatrixLike!(typeof(dftmat)));

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


    {
        auto a = dftMatrix!C(32);
        auto b = idftMatrix!C(32);
        auto c = a * b;
        auto d = b * a;
        static assert(is(typeof(c) == typeof(identity!C(32))));
        static assert(is(typeof(d) == typeof(identity!C(32))));
    }
}


struct DFTMatrixMulMM(C, Flag!"isFwd" isFwd, Flag!"applyFromLHS" applyFromLHS = Yes.applyFromLHS, Mat)
if(isMatrixLike!Mat)
{
    private alias StupidGEMM = MatrixMatrixMulGEMM!(C, DFTMatrix!(C, isFwd), Mat, ZeroMatrix!C);
    alias ElementType = StupidGEMM.ElementType;

    import std.math : log2;
    enum float exprTreeCost = Mat.exprTreeCost + EXPR_COST_N^^2 * log2(EXPR_COST_N);


    size_t length(size_t dim)() const
    {
        static if(dim == 0)
            return applyFromLHS ? _dftM.length!0 : _target.length!0;
        else
            return applyFromLHS ? _target.length!1 : _dftM.length!1;
    }


    ElementType opIndex(size_t i, size_t j) const
    {
        ElementType sum = ElementType(0);
        foreach(n; 0 .. _dftM.length!1) {
            static if(applyFromLHS)
                sum += _dftM[i, n] * _target[n, j];
            else
                sum += _target[i, n] * _dftM[n, j];
        }
        
        return sum;
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc)
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        alias F = typeof(ElementType.init.re);
        auto fftw = _dftM._fftw;

        import dffdd.math.math : fast_sqrt;
        F norm = fast_sqrt!F(_dftM.length!0);

        _target.evalTo(dst, alloc);
        foreach(j; 0 .. this.length!1) {
            static if(applyFromLHS)
                fftw.inputs!F.sliced[] = dst[0 .. $, j];
            else
                fftw.inputs!F.sliced[] = dst[j, 0 .. $];

            static if(isFwd)
                fftw.fft!F();
            else
                fftw.ifft!F();

            static if(applyFromLHS)
                dst[0 .. $, j] = fftw.outputs!F.sliced[];
            else
                dst[j, 0 .. $] = fftw.outputs!F.sliced[];
        }

        static if(isFwd)
            dst[] /= norm;
        else
            dst[] *= norm;
    }


    auto _opB_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(mat.length!0 == this.length!1)
    {
        static if(applyFromLHS)
        {
            return _dftM * (_target * mat);
        }
        else
        {
            return _target * (_dftM * mat);
        }
    }


    auto _opBR_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(mat.length!1 == this.length!0)
    {
        static if(applyFromLHS)
            return (mat * _dftM) * _target;
        else
            return (mat * _target) * _dftM;
    }


    auto _opB_Impl_(string op : "*", V)(V vec)
    if(isVectorLike!V)
    in(vec.length == this.length!1)
    {
        static if(applyFromLHS)
            return _dftM * (_target * vec);
        else
            return _target * (_dftM * vec);
    }


    auto T()()
    {
        static if(applyFromLHS)
            return _target.T * _dftM.T;
        else
            return _dftM.T * _target.T;
    }


    auto H()()
    {
        static if(applyFromLHS)
            return _target.H * _dftM.H;
        else
            return _dftM.H * _target.H;
    }


    import dffdd.math.exprtemplate;
    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*S"]));


  private:
    DFTMatrix!(C, isFwd) _dftM;
    Mat _target;
}


unittest
{
    import std.math : isClose;

    alias C = Complex!double;

    auto dftmat = DFTMatrixMulMM!(C, Yes.isFwd, Yes.applyFromLHS, Identity!C)(dftMatrix!C(4), identity!C(4)) * 2;
    auto idftmat = DFTMatrixMulMM!(C, No.isFwd, Yes.applyFromLHS, Identity!C)(idftMatrix!C(4), identity!C(4)) * 2;

    static assert(isMatrixLike!(typeof(dftmat)));

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

    auto s = slice!C(4, 4);
    dftmat.evalTo(s, theAllocator);

    assert(isClose(s[0, 0].re, 1, REPS, AEPS) && isClose(s[0, 0].im, 0, REPS, AEPS));
    assert(isClose(s[0, 1].re, 1, REPS, AEPS) && isClose(s[0, 1].im, 0, REPS, AEPS));
    assert(isClose(s[0, 2].re, 1, REPS, AEPS) && isClose(s[0, 2].im, 0, REPS, AEPS));
    assert(isClose(s[0, 3].re, 1, REPS, AEPS) && isClose(s[0, 3].im, 0, REPS, AEPS));

    assert(isClose(s[1, 0].re, 1, REPS, AEPS) && isClose(s[1, 0].im, 0, REPS, AEPS));
    assert(isClose(s[1, 1].re, 0, REPS, AEPS) && isClose(s[1, 1].im, -1, REPS, AEPS));
    assert(isClose(s[1, 2].re, -1, REPS, AEPS) && isClose(s[1, 2].im, 0, REPS, AEPS));
    assert(isClose(s[1, 3].re, 0, REPS, AEPS) && isClose(s[1, 3].im, 1, REPS, AEPS));

    assert(isClose(s[2, 0].re, 1, REPS, AEPS) && isClose(s[2, 0].im, 0, REPS, AEPS));
    assert(isClose(s[2, 1].re, -1, REPS, AEPS) && isClose(s[2, 1].im, 0, REPS, AEPS));
    assert(isClose(s[2, 2].re, 1, REPS, AEPS) && isClose(s[2, 2].im, 0, REPS, AEPS));
    assert(isClose(s[2, 3].re, -1, REPS, AEPS) && isClose(s[2, 3].im, 0, REPS, AEPS));

    assert(isClose(s[3, 0].re, 1, REPS, AEPS) && isClose(s[3, 0].im, 0, REPS, AEPS));
    assert(isClose(s[3, 1].re, 0, REPS, AEPS) && isClose(s[3, 1].im, 1, REPS, AEPS));
    assert(isClose(s[3, 2].re, -1, REPS, AEPS) && isClose(s[3, 2].im, 0, REPS, AEPS));
    assert(isClose(s[3, 3].re, 0, REPS, AEPS) && isClose(s[3, 3].im, -1, REPS, AEPS));

    idftmat.evalTo(s, theAllocator);

    assert(isClose(s[0, 0].re, 1, REPS, AEPS) && isClose(s[0, 0].im, 0, REPS, AEPS));
    assert(isClose(s[0, 1].re, 1, REPS, AEPS) && isClose(s[0, 1].im, 0, REPS, AEPS));
    assert(isClose(s[0, 2].re, 1, REPS, AEPS) && isClose(s[0, 2].im, 0, REPS, AEPS));
    assert(isClose(s[0, 3].re, 1, REPS, AEPS) && isClose(s[0, 3].im, 0, REPS, AEPS));

    assert(isClose(s[1, 0].re, 1, REPS, AEPS) && isClose(s[1, 0].im, 0, REPS, AEPS));
    assert(isClose(s[1, 1].re, 0, REPS, AEPS) && isClose(s[1, 1].im, 1, REPS, AEPS));
    assert(isClose(s[1, 2].re, -1, REPS, AEPS) && isClose(s[1, 2].im, 0, REPS, AEPS));
    assert(isClose(s[1, 3].re, 0, REPS, AEPS) && isClose(s[1, 3].im, -1, REPS, AEPS));

    assert(isClose(s[2, 0].re, 1, REPS, AEPS) && isClose(s[2, 0].im, 0, REPS, AEPS));
    assert(isClose(s[2, 1].re, -1, REPS, AEPS) && isClose(s[2, 1].im, 0, REPS, AEPS));
    assert(isClose(s[2, 2].re, 1, REPS, AEPS) && isClose(s[2, 2].im, 0, REPS, AEPS));
    assert(isClose(s[2, 3].re, -1, REPS, AEPS) && isClose(s[2, 3].im, 0, REPS, AEPS));

    assert(isClose(s[3, 0].re, 1, REPS, AEPS) && isClose(s[3, 0].im, 0, REPS, AEPS));
    assert(isClose(s[3, 1].re, 0, REPS, AEPS) && isClose(s[3, 1].im, -1, REPS, AEPS));
    assert(isClose(s[3, 2].re, -1, REPS, AEPS) && isClose(s[3, 2].im, 0, REPS, AEPS));
    assert(isClose(s[3, 3].re, 0, REPS, AEPS) && isClose(s[3, 3].im, 1, REPS, AEPS));


    auto rhs = matrix!C(4, 4, C(0));
    rhs[0, 0] = C(1, 1); rhs[0, 1] = C(2, 1); rhs[0, 2] = C(1, 3); rhs[0, 3] = C(4, 1);
    rhs[1, 0] = C(1, 2); rhs[1, 1] = C(3, 1); rhs[1, 2] = C(3, 3); rhs[1, 3] = C(4, 1);
    rhs[2, 0] = C(1, 5); rhs[2, 1] = C(4, 1); rhs[2, 2] = C(1, 3); rhs[2, 3] = C(4, 3);
    rhs[3, 0] = C(1, 1); rhs[3, 1] = C(2, 0); rhs[3, 2] = C(0, 3); rhs[3, 3] = C(4, 1);

    auto mul1 = DFTMatrixMulMM!(C, Yes.isFwd, Yes.applyFromLHS, typeof(rhs))(dftMatrix!C(4), rhs);

    import dffdd.math.exprtemplate : matrixGemm;
    auto mul2 = matrixGemm(C(1), dftMatrix!C(4), rhs, C(0), zeros!C(4, 4));

    mul1.evalTo(s, theAllocator);

    foreach(i; 0 .. 4)
        foreach(j; 0 .. 4) {
            assert(isClose(mul1[i, j].re, s[i, j].re, REPS, AEPS) && isClose(mul1[i, j].im, s[i, j].im, REPS, AEPS));
            assert(isClose(mul2[i, j].re, s[i, j].re, REPS, AEPS) && isClose(mul2[i, j].im, s[i, j].im, REPS, AEPS));
        }
}


struct DFTMatrixMulMV(C, Flag!"isFwd" isFwd, Vec)
if(isVectorLike!Vec)
{
    alias ElementType = typeof(C.init * Vec.ElementType.init);

    import std.math : log2;
    enum float exprTreeCost = Vec.exprTreeCost + EXPR_COST_N * log2(EXPR_COST_N);


    size_t length() const
    {
        return _dftM.length!0;
    }


    ElementType opIndex(size_t i) const
    {
        ElementType sum = ElementType(0);
        foreach(n; 0 .. _dftM.length!1)
            sum += _dftM[i, n] * _target[n];
        
        return sum;
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc)
    in(dst.length == this.length)
    {
        alias F = typeof(ElementType.init.re);
        auto fftw = _dftM._fftw;

        import dffdd.math.math : fast_sqrt;
        F norm = fast_sqrt!F(_dftM.length!0);

        _target.evalTo(dst, alloc);
        fftw.inputs!F.sliced[] = dst;

        static if(isFwd)
            fftw.fft!F();
        else
            fftw.ifft!F();

        dst[] = fftw.outputs!F.sliced;

        static if(isFwd)
            dst[] /= norm;
        else
            dst[] *= norm;
    }


    import dffdd.math.exprtemplate;
    mixin(definitionsOfVectorOperators(["defaults", "V+V", "V*S"]));


  private:
    DFTMatrix!(C, isFwd) _dftM;
    Vec _target;
}

unittest
{
    import std.math : isClose;

    alias C = Complex!double;
    enum float REPS = 0,        // 相対誤差でのチェックを無効
               AEPS = 1e-9;     // 絶対誤差でチェック

    auto rhs = vector!C(4, C(0));
    rhs[0] = C(1, 1); rhs[1] = C(2, 1); rhs[2] = C(1, 3); rhs[3] = C(4, 1);

    auto mul1 = DFTMatrixMulMV!(C, Yes.isFwd, typeof(rhs))(dftMatrix!C(4), rhs);
    checkIsVectorLike(mul1);
    static assert(isVectorLike!(typeof(mul1)));

    import dffdd.math.exprtemplate : vectorGemv;
    auto mul2 = vectorGemv(C(1), dftMatrix!C(4), rhs, C(0), zeros!C(4));

    auto s = slice!C(4);
    mul1.evalTo(s, theAllocator);

    foreach(i; 0 .. 4) {
        assert(isClose(mul1[i].re, s[i].re, REPS, AEPS) && isClose(mul1[i].im, s[i].im, REPS, AEPS));
        assert(isClose(mul2[i].re, s[i].re, REPS, AEPS) && isClose(mul2[i].im, s[i].im, REPS, AEPS));
    }
}



struct PermutationMatrix
{
    alias ElementType = int;
    enum float exprTreeCost = 1;


    this(in size_t[] perm, Flag!"isRev" isRev = No.isRev)
    in(checkIsUnitary(perm))
    {
        _fwdperm = perm;
        auto revperm = new size_t[perm.length];
        foreach(i, e; perm)
            revperm[e] = i;

        _revperm = revperm;
        _isRev = isRev;

        // make _p2list
        size_t[] p2list = new size_t[perm.length - 1];
        foreach(i; 0 .. perm.length - 1) {
            size_t tidx = _fwdperm[i];
            while(tidx < i)
                tidx = _fwdperm[tidx];
            
            p2list[i] = tidx;
        }
        _p2list = p2list;
    }


    size_t length(size_t dim)() const
    {
        return _fwdperm.length;
    }


    ElementType opIndex(size_t i, size_t j) const
    {
        if(!_isRev)
        {
            return _fwdperm[i] == j ? ElementType(1) : ElementType(0);
        }
        else
        {
            return _revperm[i] == j ? ElementType(1) : ElementType(0);
        }
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        dst[] = T(0);

        if(!_isRev)
        {
            foreach(i, j; _fwdperm)
                dst[i, j] = T(1);
        }
        else
        {
            foreach(i, j; _revperm)
                dst[i, j] = T(1);
        }
    }


    PermutationMatrix T() const
    {
        PermutationMatrix dst = this;
        dst._isRev = !this._isRev;
        return dst;
    }


    auto _opB_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(mat.length!0 == this.length!1)
    {
        return MulToMatrix!(M, true)(mat, this);
    }


    auto _opBR_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(mat.length!1 == this.length!0)
    {
        return MulToMatrix!(M, false)(mat, this);
    }


    auto _opB_Impl_(string op : "*", V)(V vec)
    if(isVectorLike!V)
    in(vec.length == this.length!1)
    {
        return MulToVector!V(vec, this);
    }


    alias H = T;


    import dffdd.math.exprtemplate;
    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*S"]));


  private:
    const(size_t)[] _fwdperm;
    const(size_t)[] _revperm;
    const(size_t)[] _p2list;
    bool _isRev;


    static
    bool checkIsUnitary(in size_t[] perm)
    {
        auto check = new bool[perm.length];
        foreach(i, e; perm) {
            if(e >= perm.length || check[e])
                return false;

            check[e] = true;
        }

        foreach(i, e; check)
            if(!e) return false;
        
        return true;
    }


    static
    void applyRowPermutation(T, size_t dim, SliceKind kindA)(Slice!(T*, dim, kindA) dst, const(size_t)[] p2list, bool isRev)
    {
        if(dst.length!0 < 2)
            return;

        static
        void swap_by_blas(T, SliceKind kindX, SliceKind kindY)(Slice!(T*, 1, kindX) x, Slice!(T*, 1, kindY) y)
        {
            static if(isNarrowComplex!T) {{
                import cblas;
                alias F = typeof(T.init.re);
                cblas.swap(
                    cast(cblas.blasint) x.length,
                    cast(MirComplex!F*) x.iterator,
                    cast(cblas.blasint) x._stride,
                    cast(MirComplex!F*) y.iterator,
                    cast(cblas.blasint) y._stride,
                );
            }} else {
                import mir.blas : swap;
                swap(x, y);
            }
        }

        void swap_impl(size_t i, size_t j)
        {
            static if(dim == 2) {
                static if(isBlasType!T) {
                    swap_by_blas(dst[i], dst[j]);
                } else {
                    import std.algorithm : swap;
                    foreach(k; 0 .. dst.length!1)
                        swap(dst[i, k], dst[j, k]);
                }
            } else {
                import std.algorithm : swap;
                swap(dst[i], dst[j]);
            }
        }

        import dffdd.math.linalg : isBlasType;

        if(!isRev) {
            foreach(i, j; p2list) {                
                swap_impl(i, j);
            }
        } else {
            foreach_reverse(i, j; p2list) {                
                swap_impl(i, j);
            }
        }
    }


    static struct MulToMatrix(M, bool isRowPerm)
    {
        alias ElementType = M.ElementType;
        enum float exprTreeCost = M.exprTreeCost + 1;

        size_t length(size_t dim)() const
        {
            static if((isRowPerm && dim == 0) || (!isRowPerm && dim == 1))
                return _pm.length!0;
            else
                return _target.length!dim;
        }


        ElementType opIndex(size_t i, size_t j) const
        {
            static if(isRowPerm)
            {
                auto plist = _pm._isRev ? _pm._revperm : _pm._fwdperm;
                return _target[plist[i], j];
            }
            else
            {
                auto plist = !_pm._isRev ? _pm._revperm : _pm._fwdperm;
                return _target[i, plist[j]];
            }
        }


        void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc)
        in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
        {
            _target.evalTo(dst, alloc);

            static if(isRowPerm)
                applyRowPermutation(dst, _pm._p2list, _pm._isRev);
            else
                applyRowPermutation(dst.transposed, _pm._p2list, !_pm._isRev);
        }


        auto _opB_Impl_(string op : "*", M)(M mat)
        if(isMatrixLike!M)
        in(mat.length!0 == this.length!1)
        {
            static if(isRowPerm)
            {
                return _pm * (_target * mat);
            }
            else
            {
                return _target * (_pm * mat);
            }
        }


        auto _opBR_Impl_(string op : "*", M)(M mat)
        if(isMatrixLike!M)
        in(mat.length!1 == this.length!0)
        {
            static if(isRowPerm)
                return (mat * _pm) * _target;
            else
                return (mat * _target) * _pm;
        }


        auto _opB_Impl_(string op : "*", V)(V vec)
        if(isVectorLike!V)
        in(vec.length == this.length!1)
        {
            static if(isRowPerm)
                return _pm * (_target * vec);
            else
                return _target * (_pm * vec);
        }


        auto T()()
        {
            static if(isRowPerm)
                return _target.T * _pm.T;
            else
                return _pm.T * _target.T;
        }


        auto H()()
        {
            static if(isRowPerm)
                return _target.H * _pm.H;
            else
                return _pm.H * _target.H;
        }


        import dffdd.math.exprtemplate;
        mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*S"]));

      private:
        M _target;
        PermutationMatrix _pm;
    }


    static struct MulToVector(V)
    {
        alias ElementType = V.ElementType;
        enum float exprTreeCost = V.exprTreeCost + 1;

        size_t length() const
        {
            return _target.length;
        }


        ElementType opIndex(size_t i) const
        {
            auto plist = _pm._isRev ? _pm._revperm : _pm._fwdperm;
            return _target[plist[i]];
        }


        void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc)
        in(dst.length == this.length)
        {
            _target.evalTo(dst, alloc);
            applyRowPermutation(dst, _pm._p2list, _pm._isRev);
        }


        import dffdd.math.exprtemplate;
        mixin(definitionsOfVectorOperators(["defaults", "V+V", "V*S"]));


      private:
        V _target;
        PermutationMatrix _pm;
    }
}

unittest
{
    bool checkM(A, B)(A ma, B mb)
    {
        bool ret = true;
        ret = ret && ma.length!0 == mb.length!0;
        ret = ret && ma.length!1 == mb.length!1;

        foreach(i; 0 .. ma.length!0)
            foreach(j; 0 .. ma.length!1)
                ret = ret && ma[i, j] == mb[i, j];
        
        return ret;
    }

    size_t[] perm = [1, 3, 2, 0];
    auto pm1 = PermutationMatrix(perm);
    assert(pm1[0, 0] == 0); assert(pm1[0, 1] == 1); assert(pm1[0, 2] == 0); assert(pm1[0, 3] == 0);
    assert(pm1[1, 0] == 0); assert(pm1[1, 1] == 0); assert(pm1[1, 2] == 0); assert(pm1[1, 3] == 1);
    assert(pm1[2, 0] == 0); assert(pm1[2, 1] == 0); assert(pm1[2, 2] == 1); assert(pm1[2, 3] == 0);
    assert(pm1[3, 0] == 1); assert(pm1[3, 1] == 0); assert(pm1[3, 2] == 0); assert(pm1[3, 3] == 0);

    auto slice1 = slice!size_t(4, 4);
    pm1.evalTo(slice1, matvecAllocator);
    assert(checkM(pm1, slice1));

    auto pm2 = pm1.T;
    assert(checkM(pm2, slice1.transposed));
    auto slice2 = slice!size_t(4, 4);
    pm2.evalTo(slice2, matvecAllocator);
    assert(checkM(pm2, slice2));

    auto mat = matrix!int(4, 4, 0);
    mat[0, 0] = 0; mat[1, 1] = 1; mat[2, 2] = 2; mat[3, 3] = 3;
    
    auto mul1 = pm1._opB_Impl_!"*"(mat);
    static assert(isMatrixLike!(typeof(mul1)));
    assert(mul1[0, 0] == 0); assert(mul1[0, 1] == 1); assert(mul1[0, 2] == 0); assert(mul1[0, 3] == 0);
    assert(mul1[1, 0] == 0); assert(mul1[1, 1] == 0); assert(mul1[1, 2] == 0); assert(mul1[1, 3] == 3);
    assert(mul1[2, 0] == 0); assert(mul1[2, 1] == 0); assert(mul1[2, 2] == 2); assert(mul1[2, 3] == 0);
    assert(mul1[3, 0] == 0); assert(mul1[3, 1] == 0); assert(mul1[3, 2] == 0); assert(mul1[3, 3] == 0);

    mul1.evalTo(slice1, matvecAllocator);
    assert(checkM(slice1, mul1));

    auto mul2 = mat * pm1;
    static assert(isMatrixLike!(typeof(mul2)));
    assert(mul2[0, 0] == 0); assert(mul2[0, 1] == 0); assert(mul2[0, 2] == 0); assert(mul2[0, 3] == 0);
    assert(mul2[1, 0] == 0); assert(mul2[1, 1] == 0); assert(mul2[1, 2] == 0); assert(mul2[1, 3] == 1);
    assert(mul2[2, 0] == 0); assert(mul2[2, 1] == 0); assert(mul2[2, 2] == 2); assert(mul2[2, 3] == 0);
    assert(mul2[3, 0] == 3); assert(mul2[3, 1] == 0); assert(mul2[3, 2] == 0); assert(mul2[3, 3] == 0);

    mul2.evalTo(slice1, matvecAllocator);
    assert(checkM(slice1, mul2));
}

unittest
{
    size_t[] perm = [1, 3, 2, 0];
    auto pm1 = PermutationMatrix(perm);

    auto v = vector!int(4);
    v[0] = 0; v[1] = 1; v[2] = 2; v[3] = 3;

    auto w = pm1 * v;
    static assert(isVectorLike!(typeof(w)));
    assert(w[0] == 1);
    assert(w[1] == 3);
    assert(w[2] == 2);
    assert(w[3] == 0);

    auto s = slice!int(4);
    w.evalTo(s, matvecAllocator);
    assert(s[0] == 1);
    assert(s[1] == 3);
    assert(s[2] == 2);
    assert(s[3] == 0);
}


/+
unittest
{

    Matrix!(C, Contiguous) makeCirculantIIDChannelMatrix(C)(size_t N, size_t nTaps, uint seed)
    {
        auto row = new C[N];
        row[] = C(0);
        row[0] = 1;

        foreach(i; 1 .. nTaps)
            row[$-i] = i+1;

        auto mat = matrix!C(N, N, C(0));
        foreach(i; 0 .. N) {
            mat.sliced()[i, 0 .. $] = row;
            {
                auto last = row[N-1];
                foreach_reverse(j; 1 .. N)
                    row[j] = row[j-1];

                row[0] = last;
            }
        }

        return mat;
    }

    alias C = Complex!float;
    auto chMat = makeCirculantIIDChannelMatrix!C(8, 3, 0);
    auto dftMat = forceEvaluate(dftMatrix!C(8));

    auto result1 = forceEvaluate(dftMat * chMat * dftMat.H);    // 第0列目のフーリエ変換の結果が対角に並ぶ
    auto result2 = forceEvaluate(dftMat.H * chMat * dftMat);    // 標準形（第0行目のフーリエ変換の結果が対角に並ぶ）

    import std.stdio;
    writeln(result1);
    writeln(result2);

    C[] arr = new C[8];
    foreach(i; 0 .. 8)
        arr[i] = result1[i, i];
    
    import dffdd.utils.fft;
    auto fftw = makeFFTWObject!Complex(8);
    fftw.inputs!float[] = arr[];

    fftw.ifft!float();
    writeln(fftw.outputs!float);

    foreach(i; 0 .. 8)
        arr[i] = result2[i, i];

    fftw.inputs!float[] = arr[];

    fftw.ifft!float();
    writeln(fftw.outputs!float);
}
+/