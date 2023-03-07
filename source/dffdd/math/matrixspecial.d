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

    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
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

    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
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


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
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
    enum float exprTreeCost = 1;


    this(ArrayLike vec, MatVec matvec)
    {
        _diag = vec;
        _matvec = matvec;
    }


  static if(isMatrixLike!MatVec)
  {
    size_t length(size_t dim = 0)() const @property
    if(dim == 0 || dim == 1)
    {
        return _matvec.length!dim;
    }


    ElementType opIndex(size_t i, size_t j) const
    {
        return _matvec[i, j] * _diag[i];
    }
  }
  else
  {
    size_t length() const @property
    {
        return _matvec.length;
    }


    ElementType opIndex(size_t i) const
    {
        return _matvec[i] * _diag[i];
    }
  }


  static if(isMatrixLike!MatVec)
  {
    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        _matvec.evalTo(dst, alloc);
        foreach(i; 0 .. _diag.length)
            dst[i, 0 .. $] *= _diag[i];
    }
  }
  else
  {
    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 1, kindA) dst, ref Alloc alloc) const
    in(dst.length == this.length)
    {
        _matvec.evalTo(dst, alloc);
        foreach(i; 0 .. _diag.length)
            dst[i] *= _diag[i];
    }
  }


  static if(isMatrixLike!MatVec)
  {
    auto _opB_Impl_(string op : "*", V)(V vec)
    if(isVectorLike!V)
    in(vec.length == this.length!1)
    {
        auto mm = _matvec * vec;
        return MulDiagonal!(ArrayLike, typeof(mm))(_diag, mm);
    }


    auto _opB_Impl_(string op : "*", M)(M mat)
    if(isMatrixLike!M)
    in(mat.length!0 == this.length!1)
    {
        auto mm = _matvec * mat;
        return MulDiagonal!(ArrayLike, typeof(mm))(_diag, mm);
    }

    mixin(definitionsOfMatrixOperators(["defaults", "M+M", "M*S"]));
  }
  else
  {
    mixin(definitionsOfVectorOperators(["defaults", "V+V", "V*S"]));
  }


  private:
    ArrayLike _diag;
    MatVec _matvec;
}


auto diag(Iterator, SliceKind kind)(Slice!(Iterator, 1, kind) slice)
{
    auto ident = identity!(typeof(slice[0]))(slice.length);
    return MulDiagonal!(typeof(slice), typeof(ident))(slice, ident);
}


auto diag(E)(in E[] slice)
{
    auto ident = identity!(Unqual!(typeof(slice[0])))(slice.length);
    return MulDiagonal!(const(E)[], typeof(ident))(slice, ident);
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