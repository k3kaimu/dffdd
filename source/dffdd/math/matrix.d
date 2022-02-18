module dffdd.math.matrix;

import std.experimental.allocator;
import std.traits;


import mir.ndslice;
import dffdd.math.linalg;
import dffdd.math.vector;


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
    return Matrix!(T, Contiguous)(slice!T([row, col], init));
}


auto matrix(Iterator, SliceKind kind)(Slice!(Iterator, 2, kind) s)
{
    return matrixed(slice(s));
}


enum isMatrixLike(T) = is(typeof((T t, size_t i){
    T.ElementType e = t[i, i];
    size_t rowlen = t.length!0;
    size_t collen = t.length!1;

    auto dst = slice!(typeof(e))(rowlen, collen);
    auto allocator = std.experimental.allocator.theAllocator;
    t.evalTo(dst, allocator);
}));


struct MatrixedSlice(Iterator, SliceKind kind)
{
    alias ElementType = Unqual!(typeof(opIndex(0, 0)));

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
    const(Slice!(Iterator, 2, kind)) sliced() const { return _slice; }


    auto ref opIndex(size_t i, size_t j)
    {
        return _slice[i, j];
    }


    auto opIndex(size_t i, size_t j) const
    {
        return _slice[i, j];
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


    import dffdd.math.exprtemplate;
    mixin MatrixOperators!(["M*M", "M*V", "M*S", ".H", "M=M"]);


  private:
    Slice!(Iterator, 2, kind) _slice;
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
}
