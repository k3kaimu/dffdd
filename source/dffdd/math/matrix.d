module dffdd.math.matrix;

import std.experimental.allocator;
import std.traits;
import std.typecons : Tuple;

import mir.ndslice;
import dffdd.math.linalg;
import dffdd.math.vector;
import dffdd.math.exprtemplate : isExpressionTemplate;


MatrixedSlice!(Iterator, kind) matrixed(Iterator, SliceKind kind)(Slice!(Iterator, 2, kind) slice)
{
    return MatrixedSlice!(Iterator, kind)(slice);
}


Matrix!(T, Contiguous) matrix(T)(size_t row, size_t col)
{
    return matrix!T(row, col, T.init);
}


Matrix!(T, Contiguous) matrix(T)(size_t row, size_t col, T init)
{
    return MatrixedSlice!(T*, Contiguous)(slice!T([row, col], init));
}


auto matrix(Iterator, SliceKind kind)(Slice!(Iterator, 2, kind) s)
{
    auto ret = matrix!(typeof(s[0][0]))(s.length!0, s.length!1);
    ret.sliced[] = s;
    return ret;
}


auto makeMatrix(T, Alloc)(ref Alloc alloc, size_t row, size_t col)
{
    T[] arr = alloc.makeArray!(T)(row * col);
    T[] tmp = alloc.makeArray!(T)(row * col);
    return MatrixedSlice!(T*, Contiguous)(arr.sliced(row, col), tmp.sliced(row, col));
}


auto makeMatrix(T, Alloc)(ref Alloc alloc, size_t row, size_t col, T init)
{
    T[] arr = alloc.makeArray!(T)(row * col);
    T[] tmp = alloc.makeArray!(T)(row * col);
    arr[] = init;
    return MatrixedSlice!(T*, Contiguous)(arr.sliced(row, col), tmp.sliced(row, col));
}


auto disposeMatrix(T, Alloc, SliceKind kindA)(ref Alloc alloc, ref Matrix!(T*, kindA) mat)
{
    alloc.dispose(mat._slice.iterator);
    alloc.dispose(mat._tmp.iterator);
}


enum isMatrixLike(T) = isExpressionTemplate!T && is(typeof((T t, size_t i){
    T.ElementType e = t[i, i];
    size_t rowlen = t.length!0;
    size_t collen = t.length!1;

    auto dst = slice!(typeof(e))(rowlen, collen);
    auto allocator = std.experimental.allocator.theAllocator;
    t.evalTo(dst, allocator);
}));


void checkIsMatrixLike(T)(T _)
{
    auto dg = (size_t i, T t){
        T.ElementType e = t[i, i];
        size_t rowlen = t.length!0;
        size_t collen = t.length!1;

        auto dst = slice!(typeof(e))(rowlen, collen);
        auto allocator = std.experimental.allocator.theAllocator;
        t.evalTo(dst, allocator);
    };
}


struct MatrixedSlice(Iterator, SliceKind kind)
{
    alias ElementType = Unqual!(typeof(Slice!(Iterator, 2, kind).init[0, 0]));
    enum float exprTreeCost = 0;

    this(Slice!(Iterator, 2, kind) s)
    {
        _slice = s;
    }


    size_t length(size_t dim = 0)() const @property
    if(dim == 0 || dim == 1)
    {
        return _slice.length!dim;
    }


    Slice!(Iterator, 2, kind) sliced() { return _slice; }
    auto sliced() const { return _slice.lightConst; }


    auto lightConst() const
    {
        return this.sliced.lightConst.matrixed;
    }


    auto ref opIndex(size_t i, size_t j)
    {
        return _slice[i, j];
    }


    auto opIndex(size_t i, size_t j) const
    {
        return _slice[i, j];
    }


    auto rowVec(size_t r)
    {
        auto v = _slice[r, 0 .. $].vectored;
        setTMM(v);
        return v;
    }


    auto colVec(size_t c)
    {
        auto v = _slice[0 .. $, c].vectored;
        setTMM(v);
        return v;
    }


    void evalTo(T, SliceKind kindA, Alloc)(Slice!(T*, 2, kindA) dst, ref Alloc alloc) const
    in(dst.length!0 == this.length!0 && dst.length!1 == this.length!1)
    {
        dst[] = _slice.lightConst.as!(const T);
    }


    auto T()
    {
        return _slice.transposed.matrixed;
    }


    auto T() const
    {
        return _slice.lightConst.transposed.matrixed;
    }


    auto _op_partial_impl_(Tuple!(size_t, size_t) rs, Tuple!(size_t, size_t) cs)
    {
        return _slice[rs[0] .. rs[1], cs[0] .. cs[1]].matrixed;
    }


    auto _op_partial_impl_(Tuple!(size_t, size_t) rs, size_t j)
    {
        return _slice[rs[0] .. rs[1], j].vectored;
    }


    auto _op_partial_impl_(size_t i, Tuple!(size_t, size_t) cs)
    {
        return _slice[i, cs[0] .. cs[1]].vectored;
    }


    import dffdd.math.exprtemplate;
    mixin(definitionsOfMatrixOperators(["defaults", "M*M", "M+M", "M*V", "M*S", ".H"]));


    static if(is(Iterator == ElementType*))
    {
        alias noalias = noAliasCopy;


        auto noAliasCopy(V)(V mat)
        if(isMatrixLike!V)
        in(mat.length!0 == this.length!0 && mat.length!1 == this.length!1)
        {
            if(_tmm is null) _tmm = new TempMemoryManager(1024);
            auto tmms = _tmm.saveState;
            mat.evalTo(_slice, _tmm);
        }


        auto opSliceAssign(M)(M mat)
        if(isMatrixLike!M)
        in(mat.length!0 == this.length!0 && mat.length!1 == this.length!1)
        {
            if(_tmm is null) _tmm = new TempMemoryManager(1024);
            auto tmms = _tmm.saveState;
            auto tmp = _tmm.makeSlice!ElementType(_slice.shape);
            mat.evalTo(tmp, _tmm);
            _slice[] = tmp;
        }


        auto opSliceOpAssign(string op, M)(M mat)
        if(isMatrixLike!M && (op == "+" || op == "-"))
        in(mat.length!0 == this.length!0 && mat.length!1 == this.length!1)
        {
            if(_tmm is null) _tmm = new TempMemoryManager(1024);
            auto tmms = _tmm.saveState;
            auto tmp = _tmm.makeSlice!ElementType(_slice.shape);
            mat.evalTo(tmp, matvecAllocator);

            static if(op == "+")
                _slice[] += tmp;
            else
                _slice[] -= tmp;
        }


        auto opSliceAssign(S)(S scalar)
        if(!isMatrixLike!S && !isVectorLike!S && is(S : ElementType))
        {
            this._slice[] = cast(ElementType)(scalar);
        }


        auto opSliceOpAssign(string op, S)(S scalar)
        if(!isMatrixLike!S && !isVectorLike!S && (op == "+" || op == "-" || op == "*" || op == "/"))
        {
            static if(op == "+")
                this._slice[] += scalar;
            else static if(op == "-")
                this._slice[] -= scalar;
            else static if(op == "*")
                this._slice[] *= scalar;
            else static if(op == "/")
                this._slice[] /= scalar;
            else static assert(0);
        }
    }


  private:
    Slice!(Iterator, 2, kind) _slice;

    static if(is(Iterator == ElementType*))
    {
        TempMemoryManager _tmm;
    }

    void setTMM(X)(ref X m)
    {
        static if(is(Iterator == ElementType*) && is(typeof(m._tmm)))
            m._tmm = this._tmm;
    }
}


alias Matrix(T, SliceKind kind) = MatrixedSlice!(T*, kind);


unittest
{
    auto mat1 = matrix!int(2, 2, 1);
    assert(mat1.length!0 == 2);
    assert(mat1[0, 0] == 1);
    assert(mat1[0, 1] == 1);
    assert(mat1[1, 0] == 1);
    assert(mat1[1, 1] == 1);

    mat1[0, 0] = -1;
    assert(mat1[0, 0] == -1);

    auto mat2 = iota([4], 0, 1).map!"a*2".sliced(2, 2).matrix;
    assert(mat2[0, 0] == 0);
    assert(mat2[0, 1] == 2);
    assert(mat2[1, 0] == 4);
    assert(mat2[1, 1] == 6);

    auto s = slice!int(2, 2);
    mat2.evalTo(s, theAllocator);
    assert(s[0, 0] == 0);
    assert(s[0, 1] == 2);
    assert(s[1, 0] == 4);
    assert(s[1, 1] == 6);

    auto p = mat2[0 .. 1, 0 .. 1];
    assert(p[0, 0] == 0);
}


