module dffdd.utils.file;

import std.stdio;
import std.traits;


/++
ファイルから型Eの配列データを読み込みます．
this.front はbufferSize個の要素を持つバッファを返しますが，popFront()はそのうちのstride個の要素しか更新しません．
+/
struct RawDataBufferedReader(E)
if(is(Unqual!E == E))
{

    this(string filename, size_t bufferSize, size_t stride = 0)
    {
        if(stride == 0) stride = bufferSize;
        this(File(filename, "rb"), bufferSize, stride);
    }


    this(File file, size_t bufferSize, size_t stride = 0)
    {
        if(stride == 0) stride = bufferSize;
        _file = file;
        _buffer = new E[bufferSize];
        _stride = stride;

        auto readed = _file.rawRead(_buffer);
        if(readed.length != _buffer.length)
            _empty = true;
    }


    const(E)[] front() const @property
    {
        return _buffer;
    }


    void popFront()
    {
        foreach(i; _stride .. _buffer.length)
            _buffer[i - _stride] = _buffer[i];

        auto readed = _file.rawRead(_buffer[$ - _stride .. $]);
        if(readed.length != _stride)
            _empty = true;
    }


    bool empty() const @property
    {
        return _empty;
    }


  private:
    File _file;
    E[] _buffer;
    size_t _stride;
    bool _empty = false;
}


unittest 
{
    import std.file : remove, write;

    immutable filename = "remove_me.dat";
    scope(exit) remove(filename);

    alias E = int;
    E[] data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    write(filename, data);

    auto reader1 = RawDataBufferedReader!E(filename, 3, 2);
    assert(!reader1.empty);
    assert(reader1.front == [1, 2, 3]);
    reader1.popFront();
    assert(!reader1.empty);
    assert(reader1.front == [3, 4, 5]);
    reader1.popFront();
    assert(!reader1.empty);
    assert(reader1.front == [5, 6, 7]);
    reader1.popFront();
    assert(!reader1.empty);
    assert(reader1.front == [7, 8, 9]);
    reader1.popFront();
    assert(reader1.empty);


    auto reader2 = RawDataBufferedReader!E(filename, 100, 1);
    assert(reader2.empty);
}

