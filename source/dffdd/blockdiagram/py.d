module dffdd.blockdiagram.py;

import core.thread;



struct MsgpackRPCBlock(I, O, string sendName, string recvName, size_t iTypeSize = I.sizeof, size_t oTypeSize = O.sizeof)
{
    static
    auto makeBlock(R)(R r, TCPClient client)
    if(isInputRange!R && is(ElementType!R == I))
    {
        import std.typecons;

        return RefCounted!(OFDMModulatorImpl!R)(r, client);
    }


    static
    struct OFDMModulatorImpl(R)
    {
        this(R r, TCPClient client)
        {
            _r = r;
            _client = client;
            _empty = false;
            _i = 127;
            _buf.length = 128;
            popFront();
        }


        O front() const @property { return _buf[_i]; }
        bool empty() const @property { return _empty; }

        void popFront()
        {
            ++_i;

            if(_i == _buf.length){
                _i = 0;
                _buf.length = 0;
                generate();
            }
        }

      private:
        R _r;
        TCPClient _client;
        bool _empty;
        size_t _i;
        O[] _buf;

        void generate()
        {
            bool sent;
            size_t cnt;
            while(1){
                auto rep = _client.call!(O[])(recvName, 2048 / oTypeSize);
                if(rep.empty && !sent){
                    I[2048 / iTypeSize] _arr_;
                    I[] sends = _arr_[];

                    foreach(i; 0 .. _arr_.length+1){
                        if(_r.empty || i == _arr_.length+1){
                            sends = sends[0 .. i];
                            break;
                        }

                        sends[i] = _r.front;
                        _r.popFront();
                    }

                    if(sends.length == 0){
                        assert(_r.empty);
                        ++cnt;
                        if(cnt == 1024){
                            _empty = true;
                            return;
                        }
                        Thread.sleep(dur!"msecs"(10));
                        continue;
                    }
                    _client.notify(sendName, sends);
                }
                else if(rep.empty && sent){
                    break;
                }else{
                    sent = true;
                    _buf ~= rep;
                }
            }
        }
    }
}


//final class PythonEnvironment
//{
//    void addImport(string code)
//    {
//        _imports[code] = true;
//    }


//    void addDef(string tag, string code)
//    {
//        _defs ~= code;
//    }


//    void addStmts(string code)
//    {
//        _mainstmts ~= code;
//    }


//    bool[string] _imports;
//    string[string] _defs;
//    string[] _mainstmts;
//}


//auto defGenericPushBlock(PythonEnvironment pe)
//{
//    pe.addImport("from gnuradio import gr");
//    pe.addDef("generic_push_block", q{
//class generic_push_block(gr.sync_block):
//    def __init__(self, tp, fn):
//        self.fn = fn
//        gr.sync_block.__init__(
//            self,
//            name = "genric_push_block",
//            in_sig = [tp],
//            out_sig = []
//        )

//    def work(self, input_items, output_items):
//        self.fn(input_items[0])
//        return len(input_items[0])
//});
//}


//auto defGenericPullBlock(PythonEnvironment pe)
//{
//    pe.addImport("from gnuradio import gr");
//    pe.addDef("generic_pull_block", q{
//class generic_pull_block(gr.sync_block):
//    def __init__(self, tp, fn):
//        self.fn = fn
//        gr.sync_block.__init__(
//            self,
//            name = "generic_pull_block",
//            in_sig = [],
//            out_sig = [tp]
//        )

//    def work(self, input_items, output_items):
//        arr = self.fn()
//        l = min(len(arr), len(output_items[0]))
//        output_items[0][0:l] = arr[0:l]
//        return l
//});


//    static struct MakeGenericPullBlock
//    {

//    }
//}


//auto defOFDMMod(PythonEnvironment pe)

//import std.algorithm;
//import std.range;
//import std.string;
//import std.traits;
//import std.conv;

//import pyd.pyd,
//       pyd.embedded;



//class PythonMsgpackRPCHelper
//{
//    this() { }
//}


//class PythonGNURadioMsgpackRPCServer : PythonMsgpackRPCHelper
//{
//    this(string src) {}
//}


//class MsgpackRPCHelper(Interfaces...)
//{
//    void connect(string host, ushort port)
//    {

//    }


//  private:
//    Interfaces _interfaces;
//}

//private
//auto makeProcessManager(Pid pid)
//{
//    static struct Proceess
//    {
//        this(Pid pid)
//        {
//            this._pid = pid;
//        }


//        ~this()
//        {
//            kill(_pid);
//            wait(_pid);
//        }
//    }
//}


//auto make





/+
shared static this()
{
    py_init();
}


template numpyType(T)
{
  static if(isIntegral!T && isUnsigned!T)
    enum string numpyType = "uint" ~ to!string(T.sizeof*8);
  else static if(isIntegral!T && isSigned!T)
    enum string numpyType = "int" ~ to!string(T.sizeof*8);
  else static if(isFloatingPoint!T)
    enum string numpyType = "float" ~ to!string(T.sizeof*8);
  else static if(is(typeof((T e){ auto r = e.re, i = e.im; })))
    enum string numpyType = "complex" ~ to!string(T.sizeof*8);
  else
    static assert(0);
}


final class PythonGRBlock(R, E)
{
    this(R r, string code, size_t inputBufSize = 1024)
    {
        import std.range;

        _r = r;
        _inpBufSize = inputBufSize;
        _py = new InterpContext();

        import std.stdio;
        writeln(code.replace("%(inputType)", numpyType!(ElementType!R))
                                 .replace("%(outputType)", numpyType!(E))
                                 );
        _py.py_stmts(code.replace("%(inputType)", numpyType!(ElementType!R))
                                 .replace("%(outputType)", numpyType!(E)));
        
        //_py.ioobj = _py.py_eval("IOObject()", "ioobj");
        //_py.topBlock = _py.py_eval("TopBlock(ioobj)", "topBlock");
        //_py.py_stmts("topBlock.start()");
        //_py.topBlock.start();
    }


    E front() const @property { return _buf[_i]; }
    bool empty() const @property { return _empty; }


    void popFront()
    {
        ++_i;
        if(_i == _buf.length)
            generate();
    }


  private:
    R _r;
    size_t _inpBufSize;
    InterpContext _py;
    bool _empty;
    size_t _i;
    E[] _buf;


    void generate()
    {
        void check()
        {
            foreach(e; _py.ioobj.popOuputs()){
                _buf ~= e.to_d!E;
            }
        }

        _i = 0;
        _buf.length = 0;
        check();

        if(_i != 0) return;

        if(_r.empty){
            foreach(i; 0 .. 1000){
                Thread.sleep(dur!"msecs"(1));
                check();
                if(_buf.length)
                    return;
            }
        }else{
            // Pythonに次の入力を渡す
            ElementType!R[] buf;
            foreach(i; 0 .. _inpBufSize){
                if(!_r.empty){
                    buf ~= _r.front;
                    _r.popFront();
                }else
                    break;
            }
            _py.ioobj.pushInputs(buf);
        }
    }
}


private:

immutable string[] pyImports = [
    "import time",
    "from threading import Lock",
    "import numpy",
    "from gnuradio import gr",
    "from grc_gnuradio import wxgui as grc_wxgui",
];

// params: inputType, outputType, blockDefs, blockConnections
immutable string topBlockCode = `
class TopBlock(grc_wxgui.top_block_gui):
    def __init__(self, ioobj):
        grc_wxgui.top_block_gui.__init__(self, title="block")

        self.ioobj = ioobj
        self.src_block = generic_pull_block(numpy.%(inputType), lambda: getInputs(self.ioobj))
        self.dst_block = generic_push_block(numpy.%(outputType), lambda arr: self.ioobj.pushOutputs(arr))

%(blockDefs)
%(blockConnections)
`;


immutable string ioObjectCode = `
class generic_push_block(gr.sync_block):
    def __init__(self, tp, fn):
        self.fn = fn
        gr.sync_block.__init__(
            self,
            name = "genric_push_block",
            in_sig = [tp],
            out_sig = []
        )

    def work(self, input_items, output_items):
        self.fn(input_items[0])
        return len(input_items[0])


class generic_pull_block(gr.sync_block):
    def __init__(self, tp, fn):
        self.fn = fn
        gr.sync_block.__init__(
            self,
            name = "generic_pull_block",
            in_sig = [],
            out_sig = [tp]
        )

    def work(self, input_items, output_items):
        arr = self.fn()
        l = len(arr)
        output_items[0][0:l] = arr
        return l


def getInputs(ioobj):
    for i in xrange(10000):
        arr = ioobj.popInputs()
        if(len(arr) != 0):
            return arr
        time.sleep(0.0001)

class IOObject:
    def __init__(self):
        import threading
        self.lock = threading.Lock()
        self.inputs = []
        self.outputs = []

    def popInputs(self):
        with self.lock:
            ret = self.inputs
            self.inputs = []
            return ret

    def popOutputs(self):
        with self.lock:
            ret = self.outputs
            self.outputs = []
            return ret

    def pushInputs(self, arr):
        with self.lock:
            self.inputs.append(arr)

    def pushOutputs(self, arr):
        with self.lock:
            self.outputs.append(arr)
`;


string removeIndent(string lines, string ln = "\n")
{
    import std.stdio;

    //writeln(lines);
    auto ls = lines.splitLines().map!chomp.filter!"a.length".array();
    //writeln(ls);
    auto ids = ls.map!(a => a.countUntil!"a != ' '").fold!min(size_t.max);
    if(ids == size_t.max) return lines;

    foreach(i, ref e; ls){
        e = e[ids .. $];
        writefln("%s: %s",i+1, e);
    }

    return ls.join(ln);
}


string addIndent(string lines, size_t level, string ln = "\n")
{
    auto ls = lines.splitLines.map!chomp;

    string dst;
    foreach(e; ls){
        foreach(i; 0 .. level)
            dst ~= ' ';
        dst ~= e;
        dst ~= ln;
    }

    return dst;
}
+/