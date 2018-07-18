module dffdd.blockdiagram.block;

import dffdd.blockdiagram.utils : isConverter;

import std.algorithm;
import std.range;
import std.traits;

interface IOutputTerminal(T)
{
    T value() @property;
    ComputeBlock parent() @property;
}


abstract class ComputeBlock
{
    this() {}


    void resetComputed() @property
    {
        _isComputed = false;
        foreach(e; this.dependencies)
            e.resetComputed();
    }


    bool isComputed() const @property
    {
        return _isComputed;
    }


    void compute()
    {
        if(!this.isComputed){
            foreach(e; this.dependencies)
                e.compute();

            computeImpl();
        }

        _isComputed = true;
    }


    void addDependency(ComputeBlock block)
    {
        _dependencies ~= block;
    }


    auto dependencies()
    {
        static
        ComputeBlock torvalue(ComputeBlock c) { return c; }

        return _dependencies.map!torvalue;
    }


  protected:
    void computeImpl();


  private:
    bool _isComputed;
    ComputeBlock[] _dependencies;
}


unittest
{
    int a;

    class TestBlock : ComputeBlock
    {
        this(ComputeBlock[] deps)
        {
            super();
            foreach(d; deps)
                addDependency(d);
        }

        override
        void computeImpl()
        {
            ++a;
        }
    }

    auto block1 = new TestBlock([]);
    auto block2 = new TestBlock([block1]);

    assert(!block1.isComputed);
    assert(!block2.isComputed);

    assert(a == 0);
    block2.compute();
    assert(a == 2);
    assert(block1.isComputed);
    assert(block2.isComputed);
    block2.compute();
    assert(a == 2);
    block2.resetComputed();
    assert(!block1.isComputed);
    assert(!block2.isComputed);
}



class BufferedOutputTerminal(E) : IOutputTerminal!(const(E)[])
{
    this(BufferedComputeBlock!E block, size_t idx)
    {
        _block = block;
        _idx = idx;
    }

    override
    const(E)[] value() const @property
    {
        return _block.outputBuffer(_idx);
    }

    override
    ComputeBlock parent() @property
    {
        return _block;
    }


    BufferedComputeBlock!E _block;
    size_t _idx;
}


abstract class BufferedComputeBlock(E) : ComputeBlock
{
    this(size_t outputDim)
    {
        _numOfOutputs = outputDim;
    }


    BufferedOutputTerminal!E output(size_t idx) @property
    {
        return new BufferedOutputTerminal!E(this, idx);
    }


    BufferedOutputTerminal!E output() @property
    {
        return this.output(0);
    }


    size_t numOfOutputs() const @property { return _numOfOutputs; }


  protected:
    const(E)[] outputBuffer(size_t) const;

  private:
    size_t _numOfOutputs;
}


unittest
{
    int[] buffer = new int[10];

    class TestBlock1 : BufferedComputeBlock!int
    {
        this()
        {
            super(1);
        }


        override
        void computeImpl() {
            ++buffer[0];
        }

        override
        const(int)[] outputBuffer(size_t idx) const { return buffer; }
    }

    auto block = new TestBlock1();
    assert(buffer[0] == 0);
    assert(block.numOfOutputs == 1);
    assert(!block.isComputed);
    assert(block.output.value is buffer);
    assert(block.output.parent is block);

    block.compute();
    assert(buffer[0] == 1);
    assert(block.isComputed);
}


class InfiniteRangeSource(R, E = ElementType!R) : BufferedComputeBlock!E
if(isInputRange!R && isInfinite!R && is(ElementType!R : E))
{
    this(size_t len, R range)
    {
        super(1);
        _outputs.length = len;
        _range = range;
    }


    override
    void computeImpl()
    {
        foreach(ref e; _outputs){
            e = _range.front;
            _range.popFront();
        }
    }


    override
    const(E)[] outputBuffer(size_t idx) const
    {
        return _outputs;
    }


  private:
      E[] _outputs;
      R _range;
}


auto makeInfiniteRangeSource(R)(R range, size_t buflen)
{
    return new InfiniteRangeSource!R(buflen, range);
}

unittest
{
    import std.range;
    auto source = makeInfiniteRangeSource(sequence!"n", 10);
    auto output = source.output;

    assert(!source.isComputed);
    assert(output.value.length == 10);
    assert(output.parent is source);
    assert(source.dependencies.length == 0);

    source.compute();
    assert(output.value == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    assert(source.isComputed);

    source.resetComputed();
    assert(!source.isComputed);
    source.compute();
    assert(output.value == [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]);
}


class Adder(E) : BufferedComputeBlock!E
{
    this(size_t len, BufferedOutputTerminal!E[] sources...)
    {
        super(1);
        _sources = sources.dup;
        foreach(e; sources)
            this.addDependency(e.parent);
    
        _outputs.length = len;
    }


    override
    void computeImpl()
    {
        _outputs[] = E(0);

        foreach(terminal; _sources) {
            foreach(i, e; terminal.value)
                _outputs[i] += e;
        }
    }


    override
    const(E)[] outputBuffer(size_t idx) const
    {
        return _outputs;
    }


  private:
    E[] _outputs;
    BufferedOutputTerminal!E[] _sources;
    ComputeBlock[] _dependencies;
}

unittest
{
    auto src1 = sequence!"n".makeInfiniteRangeSource(5);
    auto src2 = sequence!"n+1".makeInfiniteRangeSource(5);

    auto adder = new Adder!ulong(5, src1.output, src2.output);
    adder.compute();
    assert(src1.isComputed);
    assert(src2.isComputed);
    assert(adder.isComputed);
    assert(adder.output.value == [1, 3, 5, 7, 9]);

    adder.resetComputed();
    assert(!src1.isComputed);
    assert(!src2.isComputed);
    assert(!adder.isComputed);
    adder.compute();
    assert(adder.output.value == [11, 13, 15, 17, 19]);
    assert(src1.isComputed);
    assert(src2.isComputed);
    assert(adder.isComputed);
}


class ConverterBlock(Conv) : BufferedComputeBlock!(Conv.OutputElementType)
{
    alias IType = Conv.InputElementType;
    alias OType = Conv.OutputElementType;


    this(size_t buflen, BufferedOutputTerminal!IType src, Conv conv)
    {
        super(1);
        this.addDependency(src.parent);
        _src = src;
        _outputs.length = buflen;
        _conv = conv;
    }


    override
    void computeImpl()
    {
        _conv(_src.value, _outputs);
    }


    override
    const(OType)[] outputBuffer(size_t idx) const
    {
        return _outputs;
    }


    ref Conv converter() @property
    {
        return _conv;
    }


  private:
    OType[] _outputs;
    Conv _conv;
    BufferedOutputTerminal!IType _src;
}


unittest
{
    static struct MultiplyConverter
    {
        alias InputElementType = int;
        alias OutputElementType = long;

        int _a;
        this(int a) { _a = a; }
        void opCall(in InputElementType[] input, OutputElementType[] output)
        {
            foreach(i; 0 .. input.length)
                output[i] = input[i] * _a;
        }
    }

    auto mc = MultiplyConverter(10);
    auto src = sequence!"n".map!"cast(int)a".makeInfiniteRangeSource(3);
    auto applied = new ConverterBlock!(MultiplyConverter)(3, src.output, mc);
    applied.compute();
    assert(applied.output.value == [0, 10, 20]);
}



class ModulatorBlockImpl(Mod, bool isMod)
    : BufferedComputeBlock!(Select!(isMod, Mod.OutputElementType, Mod.InputElementType))
{
  static if(isMod)
  {
    alias IType = Mod.InputElementType;
    alias OType = Mod.OutputElementType;
  }
  else
  {
    alias OType = Mod.InputElementType;
    alias IType = Mod.OutputElementType;
  }


    this(size_t buflen, BufferedOutputTerminal!IType src, Mod mod)
    {
        super(1);
        this.addDependency(src.parent);
        _src = src;
        _outputs.length = buflen;
        _mod = mod;
    }


    override
    void computeImpl()
    {
        static if(isMod)
            _mod.modulate(_src.value, _outputs);
        else
            _mod.demodulate(_src.value, _outputs);
    }


    override
    const(OType)[] outputBuffer(size_t idx) const
    {
        return _outputs;
    }


    ref Mod converter() @property
    {
        return _mod;
    }


  private:
    OType[] _outputs;
    Mod _mod;
    BufferedOutputTerminal!IType _src;
}


auto makeModulatorBlock(Mod, Src)(Src src, size_t len, Mod mod)
if(is(Src == BufferedOutputTerminal!(Mod.InputElementType)))
{
    return new ModulatorBlockImpl!(Mod, true)(len, src, mod);
}


auto makeDemodulatorBlock(Mod, Src)(Src src, size_t len, Mod mod)
if(is(Src == BufferedOutputTerminal!(Mod.OutputElementType)))
{
    return new ModulatorBlockImpl!(Mod, false)(len, src, mod);
}
