module dffdd.blockdiagram.inst.file;


import std.stdio;

auto fileSink(alias func)(string filename)
{
    static struct FileSink
    {
        void put(E)(E e)
        {
            func(_file, e);
        }

      private:
        File _file;
    }


    FileSink sink;
    sink._file = File(filename, "w");
    return sink;
}


auto fileSink(Fn)(string filename, Fn fn)
{
    static struct FileSink
    {
        void put(E)(E e)
        {
            _fn(_file, e);
        }

      private:
        File _file;
        Fn _fn;
    }

    FileSink sink;
    sink._file = File(filename, "w");
    sink._fn = fn;
    return sink;
}


auto binaryFileSink(string filename)
{
    static
    void rawWriter(E)(File file, E e)
    {
        file.rawWrite((cast(ubyte*)&e)[0 .. E.sizeof]);
    }


    return fileSink!rawWriter(filename);
}



auto formattedFileSink(string filename, string format)
{
    static struct Fn
    {
        void opCall(E)(File file, E e)
        {
            file.writef(_fmt, e);
        }

      private:
        string _fmt;
    }

    Fn fn;
    fn._fmt = format;

    return fileSink(filename, fn);
}


auto binaryFileReader(C)(string filename, size_t offset = 0)
{
    static struct Result
    {
        C front() const @property { return _buf[_i]; }
        bool empty() const @property { return _empty; }


        void popFront()
        {
            ++_i;
            if(_i >= _buf.length){
                _i = 0;
                rebuffering();
            }
        }


      private:
        File _file;
        C[] _buf;
        bool _empty;
        size_t _i;


        void rebuffering()
        {
            _buf = _file.rawRead(_buf);
            if(_buf.length == 0)
                _empty = true;
        }
    }


    auto res = Result(File(filename, "r"), new C[1024*1024], false, 0);
    if(offset) res._file.seek(C.sizeof * offset);
    res.rebuffering();
    return res;
}
