module dffdd.math.lapack;

import std.algorithm : min, max;
import mir.lapack;
import lapack;
import mir.ndslice;

import dffdd.math.complex;

///
size_t gebrd(T)(
    Slice!(T*, 2, Canonical) a,
    Slice!(RealPartType!T*, 1, Contiguous) d,
    Slice!(RealPartType!T*, 1, Contiguous) e,
    Slice!(T*, 1, Contiguous) tauq,
    Slice!(T*, 1, Contiguous) taup,
    Slice!(T*, 1, Contiguous) work,
    )
in
{
    assert(d.length == min(a.length!0, a.length!1));
    assert(e.length == min(a.length!0, a.length!1) - 1);
    assert(tauq.length == min(a.length!0, a.length!1));
    assert(taup.length == min(a.length!0, a.length!1));
    assert(work.length >= max(1, a.length!0, a.length!1));
}
do
{
    lapackint m = cast(lapackint) a.length!1;
    lapackint n = cast(lapackint) a.length!0;
    lapackint lda = cast(lapackint) a._stride.max(1);
    lapackint lwork = cast(lapackint) work.length;
    lapackint info;

    lapack.gebrd_(m, n, a.iterator, lda, d.iterator, e.iterator, tauq.iterator, taup.iterator, work.iterator, lwork, info);

    assert(info >= 0);
    return info;
}



size_t orgbr(T)(
    char vect,
    size_t k,
    Slice!(T*, 2, Canonical) a,
    Slice!(T*, 1, Contiguous) tau,
    Slice!(T*, 1, Contiguous) work,
)
in
{
    assert(vect == 'P' || vect == 'Q');
    assert(k >= 0);
    if(vect == 'P') {
        assert(tau.length == min(k, a.length!0));
    } else /* vect == 'Q'*/ {
        assert(tau.length == min(k, a.length!1));
    }
    assert(work.length >= max(1, min(a.length!0, a.length!1)));
}
do
{
    lapackint m = cast(lapackint) a.length!1;
    lapackint n = cast(lapackint) a.length!0;
    lapackint k_ = cast(lapackint) k;
    lapackint lda = cast(lapackint) a._stride.max(1);
    lapackint lwork = cast(lapackint) work.length;
    lapackint info;

    lapack.orgbr_(vect, m, n, k_, a.iterator, lda, tau.iterator, work.iterator, lwork, info);

    assert(info >= 0);
    return info;
}


size_t ungbr(T)(
    char vect,
    size_t k,
    Slice!(T*, 2, Canonical) a,
    Slice!(T*, 1, Contiguous) tau,
    Slice!(T*, 1, Contiguous) work,
)
in
{
    assert(vect == 'P' || vect == 'Q');
    assert(k >= 0);
    if(vect == 'P') {
        assert(tau.length == min(k, a.length!0));
    } else /* vect == 'Q'*/ {
        assert(tau.length == min(k, a.length!1));
    }
    assert(work.length >= max(1, min(a.length!0, a.length!1)));
}
do
{
    lapackint m = cast(lapackint) a.length!1;
    lapackint n = cast(lapackint) a.length!0;
    lapackint k_ = cast(lapackint) k;
    lapackint lda = cast(lapackint) a._stride.max(1);
    lapackint lwork = cast(lapackint) work.length;
    lapackint info;

    lapack.ungbr_(vect, m, n, k_, a.iterator, lda, tau.iterator, work.iterator, lwork, info);

    assert(info >= 0);
    return info;
}



///
size_t hetrd(T)(
    char uplo,
    Slice!(T*, 2, Canonical) a,
    Slice!(RealPartType!T*, 1, Contiguous) d,
    Slice!(RealPartType!T*, 1, Contiguous) e,
    Slice!(T*, 1, Contiguous) tauq,
    Slice!(T*, 1, Contiguous) work,
    )
in
{
    assert(a.length!0 == a.length!1);
    assert(d.length == a.length!0);
    assert(e.length == a.length!0 - 1);
    assert(tauq.length == a.length!0 - 1);
    assert(work.length >= max(1, a.length!0));
}
do
{
    lapackint n = cast(lapackint) a.length!0;
    lapackint lda = cast(lapackint) a._stride.max(1);
    lapackint lwork = cast(lapackint) work.length;
    lapackint info;

    lapack.hetrd_(uplo, n, a.iterator, lda, d.iterator, e.iterator, tauq.iterator, work.iterator, lwork, info);

    assert(info >= 0);
    return info;
}


///
size_t sytrd(T)(
    char uplo,
    Slice!(T*, 2, Canonical) a,
    Slice!(T*, 1, Contiguous) d,
    Slice!(T*, 1, Contiguous) e,
    Slice!(T*, 1, Contiguous) tauq,
    Slice!(T*, 1, Contiguous) work,
    )
in
{
    assert(a.length!0 == a.length!1);
    assert(d.length == a.length!0);
    assert(e.length == a.length!0 - 1);
    assert(tauq.length == a.length!0 - 1);
    assert(work.length >= max(1, a.length!0));
}
do
{
    lapackint n = cast(lapackint) a.length!0;
    lapackint lda = cast(lapackint) a._stride.max(1);
    lapackint lwork = cast(lapackint) work.length;
    lapackint info;

    lapack.sytrd_(uplo, n, a.iterator, lda, d.iterator, e.iterator, tauq.iterator, work.iterator, lwork, info);

    assert(info >= 0);
    return info;
}


size_t orgtr(T)(
    char uplo,
    Slice!(T*, 2, Canonical) a,
    Slice!(T*, 1, Contiguous) tau,
    Slice!(T*, 1, Contiguous) work,
)
in
{
    assert(uplo == 'U' || uplo == 'L');
    assert(a.length!0 == a.length!1);
    assert(tau.length == a.length!0 - 1);
    assert(work.length >= max(1, min(a.length!0, a.length!1)));
}
do
{
    lapackint n = cast(lapackint) a.length!0;
    lapackint lda = cast(lapackint) a._stride.max(1);
    lapackint lwork = cast(lapackint) work.length;
    lapackint info;

    lapack.orgtr_(uplo, n, a.iterator, lda, tau.iterator, work.iterator, lwork, info);

    assert(info >= 0);
    return info;
}


size_t ungtr(T)(
    char uplo,
    Slice!(T*, 2, Canonical) a,
    Slice!(T*, 1, Contiguous) tau,
    Slice!(T*, 1, Contiguous) work,
)
in
{
    assert(uplo == 'U' || uplo == 'L');
    assert(a.length!0 == a.length!1);
    assert(tau.length == a.length!0 - 1);
    assert(work.length >= max(1, min(a.length!0, a.length!1)));
}
do
{
    lapackint n = cast(lapackint) a.length!0;
    lapackint lda = cast(lapackint) a._stride.max(1);
    lapackint lwork = cast(lapackint) work.length;
    lapackint info;

    lapack.ungtr_(uplo, n, a.iterator, lda, tau.iterator, work.iterator, lwork, info);

    assert(info >= 0);
    return info;
}