module dffdd.blockdiagram.snitch;

import std.range;


struct Snitcher(R, O)
{
    this(R r, O o)
    {
        _r = r;
        _o = o;
        .put(_o, _r.front);
    }


    auto ref front() @property
    {
        return _r.front;
    }


    void popFront()
    {
        _r.popFront();
        if(!_r.empty) .put(_o, _r.front);
    }


    bool empty() @property
    {
        return _r.empty;
    }


  private:
    R _r;
    O _o;
}


private
struct SnitchToPointer(E)
{
    void put(X)(X a)
    {
        *_ptr = a;
    }

    private:
    E* _ptr;
}


auto connectToSnitch(R, E)(R range, E* ptr)
{
    return Snitcher!(R, SnitchToPointer!E)(range, SnitchToPointer!E(ptr));
}


unittest
{
    import std.range;

    auto r1 = iota(1, 5);
    auto ptr = new long;

    auto r2 = r1.connectToSnitch(ptr);
    assert(*ptr == 1);
    assert(r2.front == 1);
    assert(!r2.empty);
    r2.popFront();
    assert(*ptr == 2);
    assert(r2.front == 2);
    assert(!r2.empty);
    r2.popFront();
    assert(*ptr == 3);
    assert(r2.front == 3);
    assert(!r2.empty);
    r2.popFront();
    assert(*ptr == 4);
    assert(r2.front == 4);
    assert(!r2.empty);
    r2.popFront();
    assert(r2.empty);
}
