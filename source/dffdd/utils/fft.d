module dffdd.utils.fft;

import std.algorithm : swap;
import std.complex;
import std.math;
import std.numeric;
import std.stdio : File;
import std.traits : Unqual, isArray;
import std.range;
import std.meta : AliasSeq;
import std.typecons : RefCounted, Tuple, tuple;

import carbon.math;
import carbon.complex;
import carbon.memory : fastPODCopy;


enum FFTLibrary
{
    phobos,
    fftw,
}


/**
FFTを計算するためのオブジェクトです．
FFTObjectは，次のように使います．

Example:
------------
// inputをFFTした結果をoutputに格納します
void doFFT(R, FftObject)(FftObject fftobj, Complex!R[] input, Complex!R[] output)
{
    auto ip = fftobj.inputs!R;
    auto op = fftobj.outputs!R;

    ip[] = input[];
    fftobj.fft();
    output[] = op[];
}


// inputをinverse-FFTした結果をoutputに格納します
void doIFFT(R, FftObject)(FftObject fftobj, Complex!R[] input, Complex!R[] output)
{
    auto ip = fftobj.inputs!R;
    auto op = fftobj.outputs!R;

    ip[] = input[];
    fftobj.ifft();
    output[] = op[];
}
------------

fftやifftの計算は，常にobj.inputs!Rが使用されます．
また，計算結果はobj.outputs!Rに格納されます．
*/
enum bool isFFTObject(FftObj) = is(typeof((FftObj obj){
    foreach(F; AliasSeq!(float, double, real)){
        auto inpArray = obj.inputs!F;
        auto outArray = obj.outputs!F;

        obj.fft!F();
        obj.ifft!F();
    }
}));


private
void fastCopy(E)(in E[] src, E[] dst)
{
    dst[] = src[];
}


private
void fastCopy(R1, R2)(R1 src, R2 dst)
if(!isArray!R1 || !isArray!R2)
{
    foreach(i; 0 .. src.length)
        dst[i] = src[i];
}


/**

*/
void fft(F, FftObj, Ri, Ro)(FftObj obj, Ri input, Ro output)
if(isInputRange!Ri && isComplex!(Unqual!(ElementType!Ri), F)
        && isInputRange!Ro && hasAssignableElements!Ro)
in{
    static if(hasLength!Ri)
        assert(input.length == obj.inputs!F.length);

    static if(hasLength!Ro)
        assert(output.length == obj.outputs!F.length);
}
do{
    fastCopy(input, obj.inputs!F);
    obj.fft!F();
    fastCopy(obj.outputs!F, output);
}

///
unittest
{
    auto fftObj = makePhobosFFTObject!Complex(2);

    foreach(F; AliasSeq!(float, double, real))
    {
        Complex!F[] inps = new Complex!F[2],
                    outs = new Complex!F[2];

        .fft!F(fftObj, [complex!F(1, 1), complex!F(1, 1)], outs);

        assert(isClose(outs[0].re, 2));
        assert(isClose(outs[0].im, 2));
        assert(isClose(outs[1].re, 0));
        assert(isClose(outs[1].re, 0));
    }
}


/**

*/
void fftFrom(F, FftObj, Ri)(FftObj obj, Ri input)
if(isInputRange!Ri && isComplex!(Unqual!(ElementType!Ri), F))
in{
    static if(hasLength!Ri)
        assert(input.length == obj.inputs!F.length);
}
do{
    fastCopy(input, obj.inputs!F);
    obj.fft!F();
}

///
unittest
{
    auto fftObj = makePhobosFFTObject!Complex(2);

    foreach(F; AliasSeq!(float, double, real))
    {
        Complex!F[] inps = new Complex!F[2];

        .fftFrom!F(fftObj, [complex!F(1, 1), complex!F(1, 1)]);

        assert(isClose(fftObj.outputs!F[0].re, 2));
        assert(isClose(fftObj.outputs!F[0].im, 2));
        assert(isClose(fftObj.outputs!F[1].re, 0));
        assert(isClose(fftObj.outputs!F[1].re, 0));
    }
}


/**

*/
void fftTo(F, FftObj, Ro)(FftObj obj, Ro output)
if(isInputRange!Ro && hasAssignableElements!Ro)
in{
    static if(hasLength!Ro)
        assert(output.length == obj.outputs!F.length);
}
do{
    obj.fft!F();
    fastCopy(obj.outputs!F, output);
}

///
unittest
{
    auto fftObj = makePhobosFFTObject!Complex(2);

    foreach(F; AliasSeq!(float, double, real))
    {
        fftObj.inputs!F[] = [complex!F(1, 1), complex!F(1, 1)];

        Complex!F[] outs = new Complex!F[2];

        .fftTo!F(fftObj, outs);

        assert(isClose(outs[0].re, 2));
        assert(isClose(outs[0].im, 2));
        assert(isClose(outs[1].re, 0));
        assert(isClose(outs[1].re, 0));
    }
}


/**

*/
void ifft(F, FftObj, Ri, Ro)(FftObj obj, Ri input, Ro output)
if(isInputRange!Ri && isComplex!(Unqual!(ElementType!Ri), F)
        && isInputRange!Ro && hasAssignableElements!Ro)
in{
    static if(hasLength!Ri)
        assert(input.length == obj.inputs!F.length);

    static if(hasLength!Ro)
        assert(output.length == obj.outputs!F.length);
}
do{
    fastCopy(input, obj.inputs!F);
    obj.ifft!F();
    fastCopy(obj.outputs!F, output);
}

///
unittest
{
    auto fftObj = makePhobosFFTObject!Complex(2);

    foreach(F; AliasSeq!(float, double, real))
    {
        Complex!F[] inps = new Complex!F[2],
                    outs = new Complex!F[2];

        .ifft!F(fftObj, [complex!F(2, 2), complex!F(0, 0)], outs);

        assert(isClose(outs[0].re, 1));
        assert(isClose(outs[0].im, 1));
        assert(isClose(outs[1].re, 1));
        assert(isClose(outs[1].re, 1));
    }
}



/**

*/
void ifftFrom(F, FftObj, Ri)(FftObj obj, Ri input)
if(isInputRange!Ri && isComplex!(Unqual!(ElementType!Ri), F))
in{
    static if(hasLength!Ri)
        assert(input.length == obj.inputs!F.length);
}
do{
    fastCopy(input, obj.inputs!F);
    obj.ifft!F();
}

///
unittest
{
    auto fftObj = makePhobosFFTObject!Complex(2);

    foreach(F; AliasSeq!(float, double, real))
    {
        Complex!F[] inps = new Complex!F[2];

        .ifftFrom!F(fftObj, [complex!F(2, 2), complex!F(0, 0)]);

        assert(isClose(fftObj.outputs!F[0].re, 1));
        assert(isClose(fftObj.outputs!F[0].im, 1));
        assert(isClose(fftObj.outputs!F[1].re, 1));
        assert(isClose(fftObj.outputs!F[1].re, 1));
    }
}


/**

*/
void ifftTo(F, FftObj, Ro)(FftObj obj, Ro output)
if(isInputRange!Ro && hasAssignableElements!Ro)
in{
    static if(hasLength!Ro)
        assert(output.length == obj.outputs!F.length);
}
do{
    obj.ifft!F();
    fastCopy(obj.outputs!F, output);
}

///
unittest
{
    auto fftObj = makePhobosFFTObject!Complex(2);

    foreach(F; AliasSeq!(float, double, real))
    {
        fftObj.inputs!F[] = [complex!F(2, 2), complex!F(0, 0)];

        Complex!F[] outs = new Complex!F[2];

        .ifftTo!F(fftObj, outs);

        assert(isClose(outs[0].re, 1));
        assert(isClose(outs[0].im, 1));
        assert(isClose(outs[1].re, 1));
        assert(isClose(outs[1].re, 1));
    }
}



private
struct FFTWObjectImpl(alias Cpx)
{
    import core.sync.mutex : Mutex;
    import dfftw3.fftw3;
    /// mutex of fftw
    __gshared static Mutex fftwMutex;

    shared static this()
    {
        if(fftwMutex is null)
            fftwMutex = new Mutex;
    }


    this(size_t fftsize)
    {
        _fftsize = fftsize;
    }


    ~this() @nogc
    {
        if(_floatFwdPlan !is null){
            forceNoGC!fftwf_destroy_plan(_floatFwdPlan);
            _floatFwdPlan = null;
        }

        if(_doubleFwdPlan !is null){
            forceNoGC!fftw_destroy_plan(_doubleFwdPlan);
            _doubleFwdPlan = null;
        }

        if(_realFwdPlan !is null){
            forceNoGC!fftwl_destroy_plan(_realFwdPlan);
            _realFwdPlan = null;
        }

        if(_floatInvPlan !is null){
            forceNoGC!fftwf_destroy_plan(_floatInvPlan);
            _floatInvPlan = null;
        }

        if(_doubleInvPlan !is null){
            forceNoGC!fftw_destroy_plan(_doubleInvPlan);
            _doubleInvPlan = null;
        }

        if(_realInvPlan !is null){
            forceNoGC!fftwl_destroy_plan(_realInvPlan);
            _realInvPlan = null;
        }

        //static assert(0, "free in/out array!!!");
        if(_float_in !is null){
            forceNoGC!fftwf_free(_float_in.ptr);
            _float_in = null;
        }

        if(_float_out !is null){
            forceNoGC!fftwf_free(_float_out.ptr);
            _float_out = null;
        }

        if(_double_in !is null){
            forceNoGC!fftw_free(_double_in.ptr);
            _double_in = null;
        }

        if(_double_out !is null){
            forceNoGC!fftw_free(_double_out.ptr);
            _double_out = null;
        }

        if(_real_in !is null){
            forceNoGC!fftwl_free(_real_in.ptr);
            _real_in = null;
        }

        if(_real_out !is null){
            forceNoGC!fftwl_free(_real_out.ptr);
            _real_out = null;
        }
    }


    void initialize(F)()
    {
        static if(is(F == float))
        {
            if(_floatFwdPlan is null){
                auto plan = makePlanAndArray!F(_fftsize, _float_in.ptr, _float_out.ptr);

                _floatFwdPlan = plan[0];
                _floatInvPlan = plan[1];
                _float_in = plan[2][0 .. _fftsize];
                _float_out = plan[3][0 .. _fftsize];
            }
        }
        else static if(is(F == double))
        {
            if(_doubleFwdPlan is null){
                auto plan = makePlanAndArray!F(_fftsize, _double_in.ptr, _double_out.ptr);

                _doubleFwdPlan = plan[0];
                _doubleInvPlan = plan[1];
                _double_in = plan[2][0 .. _fftsize];
                _double_out = plan[3][0 .. _fftsize];
            }
        }
        else static if(is(F == real))
        {
            if(_realFwdPlan is null){
                auto plan = makePlanAndArray!F(_fftsize, _real_in.ptr, _real_out.ptr);

                _realFwdPlan = plan[0];
                _realInvPlan = plan[1];
                _real_in = plan[2][0 .. _fftsize];
                _real_out = plan[3][0 .. _fftsize];
            }
        }
    }


    private
    void fftImpl(F, int direction)()
    {
        static if(is(F == float))
        {
            static if(direction == FFTW_FORWARD)
                forceNoGC!fftwf_execute(_floatFwdPlan);
            else
                forceNoGC!fftwf_execute(_floatInvPlan);
        }
        else static if(is(F == double))
        {
            static if(direction == FFTW_FORWARD)
                forceNoGC!fftw_execute(_doubleFwdPlan);
            else
                forceNoGC!fftw_execute(_doubleInvPlan);
        }
        else static if(is(F == real))
        {
            static if(direction == FFTW_FORWARD)
                forceNoGC!fftwl_execute(_realFwdPlan);
            else
                forceNoGC!fftwl_execute(_realInvPlan);
        }
        else
            static assert(0);
    }


    Cpx!F[] inputs(F)() @property
    {
        initialize!F();

        static if(is(F == float)) return _float_in;
        else static if(is(F == double)) return _double_in;
        else static if(is(F == real)) return _real_in;
    }


    Cpx!F[] outputs(F)() @property
    {
        initialize!F();

        static if(is(F == float)) return _float_out;
        else static if(is(F == double)) return _double_out;
        else static if(is(F == real)) return _real_out;
    }


    void fft(F)()
    {
        fftImpl!(F, FFTW_FORWARD);
    }


    void ifft(F)()
    {
        fftImpl!(F, FFTW_BACKWARD);

        foreach(i; 0 .. _fftsize){
            static if(is(F == float)) _float_out[i] /= _fftsize;
            else static if(is(F == double)) _double_out[i] /= _fftsize;
            else static if(is(F == real)) _real_out[i] /= _fftsize;
        }
    }


  private:
    size_t _fftsize;

    fftwf_plan _floatFwdPlan;
    fftw_plan _doubleFwdPlan;
    fftwl_plan _realFwdPlan;

    fftwf_plan _floatInvPlan;
    fftw_plan _doubleInvPlan;
    fftwl_plan _realInvPlan;

    Cpx!float[] _float_in;
    Cpx!float[] _float_out;
    Cpx!double[] _double_in;
    Cpx!double[] _double_out;
    Cpx!real[] _real_in;
    Cpx!real[] _real_out;

    static
    auto makePlanAndArray(F)(size_t n, Cpx!F* input, Cpx!F* output) @nogc
    {
        auto newInput = input;
        auto newOutput = output;

        void checkMemory()
        {
            if(newInput is null || newOutput is null)
            {
                import core.exception;
                onOutOfMemoryError();
            }
        }

        synchronized(fftwMutex)
        {
          static if(is(F == float))
          {
            if(newInput is null) newInput = cast(Cpx!F*)forceNoGC!fftwf_malloc(F.sizeof * 2 * n);
            if(newOutput is null) newOutput = cast(Cpx!F*)forceNoGC!fftwf_malloc(F.sizeof * 2 * n);
            checkMemory();
            auto planFwd = forceNoGC!fftwf_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_FORWARD, FFTW_ESTIMATE);
            auto planInv = forceNoGC!fftwf_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_BACKWARD, FFTW_ESTIMATE);
          }
          else static if(is(F == double))
          {
            if(newInput is null) newInput = cast(Cpx!F*)forceNoGC!fftw_malloc(F.sizeof * 2 * n);
            if(newOutput is null) newOutput = cast(Cpx!F*)forceNoGC!fftw_malloc(F.sizeof * 2 * n);
            checkMemory();
            auto planFwd = forceNoGC!fftw_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_FORWARD, FFTW_ESTIMATE);
            auto planInv = forceNoGC!fftw_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_BACKWARD, FFTW_ESTIMATE);
          }
          else static if(is(F == real))
          {
            if(newInput is null) newInput = cast(Cpx!F*)forceNoGC!fftwl_malloc(F.sizeof * 2 * n);
            if(newOutput is null) newOutput = cast(Cpx!F*)forceNoGC!fftwl_malloc(F.sizeof * 2 * n);
            checkMemory();
            auto planFwd = forceNoGC!fftwl_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_FORWARD, FFTW_ESTIMATE);
            auto planInv = forceNoGC!fftwl_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_BACKWARD, FFTW_ESTIMATE);
          }

            if(planFwd is null || planInv is null)
            {
                import core.exception;
                onOutOfMemoryError();
            }

            return tuple(planFwd, planInv, newInput, newOutput);
        }
    }


    import std.traits;
    static
    auto forceNoGC(alias func, Params...)(auto ref Params params) @nogc
    if(functionLinkage!func == "C")
    {
        alias R = ReturnType!func;
        alias Ps = Parameters!func;
        static
        R funcWrapper(Ps ps)
        {
            return func(ps);
        }

        alias NoGCFunc = R function(Ps) @nogc;
        return (cast(NoGCFunc)&funcWrapper)(params);
    }
}


/**
FFTWを利用したFFTObject
*/
alias FFTWObject(alias Cpx) = RefCounted!(FFTWObjectImpl!Cpx);

///
unittest
{
    static assert(isFFTObject!(FFTWObject!(Complex)));
    static assert(isFFTObject!(FFTWObject!(std_complex_t)));
    static assert(isFFTObject!(FFTWObject!(complex_t)));
}



/**
標準ライブラリのstd.numeric.Fftを利用したFFTObject
*/
final class PhobosFFTObject(alias Cpx)
{
    this(size_t size)
    {
        _size = size;
        _impl = new Fft(size);
    }


    void initialize(F)()
    {
        static if(is(F == float))
        {
            if(_float_in is null) _float_in = new Cpx!float[_size];
            if(_float_out is null) _float_out = new Cpx!float[_size];
        }
        else static if(is(F == double))
        {
            if(_double_in is null) _double_in = new Cpx!double[_size];
            if(_double_out is null) _double_out = new Cpx!double[_size];
        }
        else static if(is(F == real))
        {
            if(_real_in is null) _real_in = new Cpx!real[_size];
            if(_real_out is null) _real_out = new Cpx!real[_size];
        }
    }


    Cpx!F[] inputs(F)() @property
    {
        initialize!F();

      static if(is(F == float))
        return _float_in;
      else static if(is(F == double))
        return _double_in;
      else static if(is(F == real))
        return _real_in;
    }


    Cpx!F[] outputs(F)() @property
    {
        initialize!F();

      static if(is(F == float))
        return _float_out;
      else static if(is(F == double))
        return _double_out;
      else static if(is(F == real))
        return _real_out;
    }


    void fft(F)()
    {
      static if(is(F == float))
      {
        auto ip = cast(Complex!F*)(_float_in.ptr);
        auto op = cast(Complex!F*)(_float_out.ptr);
      }
      else static if(is(F == double))
      {
        auto ip = cast(Complex!F*)(_double_in.ptr);
        auto op = cast(Complex!F*)(_double_out.ptr);
      }
      else static if(is(F == real))
      {
        auto ip = cast(Complex!F*)(_real_in.ptr);
        auto op = cast(Complex!F*)(_real_out.ptr);
      }

      _impl.fft(ip[0 .. _size], op[0 .. _size]);
    }


    void ifft(F)()
    {
      static if(is(F == float))
      {
        auto ip = cast(Complex!F*)(_float_in.ptr);
        auto op = cast(Complex!F*)(_float_out.ptr);
      }
      else static if(is(F == double))
      {
        auto ip = cast(Complex!F*)(_double_in.ptr);
        auto op = cast(Complex!F*)(_double_out.ptr);
      }
      else static if(is(F == real))
      {
        auto ip = cast(Complex!F*)(_real_in.ptr);
        auto op = cast(Complex!F*)(_real_out.ptr);
      }

      _impl.inverseFft(ip[0 .. _size], op[0 .. _size]);
    }


  private:
    size_t _size;
    Fft _impl;
    Cpx!float[] _float_in;
    Cpx!float[] _float_out;
    Cpx!double[] _double_in;
    Cpx!double[] _double_out;
    Cpx!real[] _real_in;
    Cpx!real[] _real_out;
}

///
unittest
{
    static assert(isFFTObject!(PhobosFFTObject!(Complex)));
    static assert(isFFTObject!(PhobosFFTObject!(std_complex_t)));
    static assert(isFFTObject!(PhobosFFTObject!(complex_t)));
}


/**
高速フーリエ変換を計算するオブジェクトを作ります．
*/
auto makeFFTObject(FFTLibrary lib = FFTLibrary.phobos, alias Cpx = Complex)(size_t size)
{
  static if(lib == FFTLibrary.phobos)
    return new PhobosFFTObject!Cpx(size);
  else static if(lib == FFTLibrary.fftw)
    return FFTWObject!Cpx(size);
}


/**
標準ライブラリで提供されているstd.numeric.Fftを利用したFFTObjectを作ります．
*/
auto makePhobosFFTObject(alias Cpx = Complex)(size_t size)
{
    return makeFFTObject!(FFTLibrary.phobos, Cpx)(size);
}


/**
FFTWを利用したFFTObjectを作ります．
*/
auto makeFFTWObject(alias Cpx = Complex)(size_t size)
{
    return makeFFTObject!(FFTLibrary.fftw, Cpx)(size);
}


unittest
{
    foreach(makeObjFunc; AliasSeq!(makePhobosFFTObject, makeFFTWObject))
        foreach(F; AliasSeq!(float, double, real))
        {
            auto fft = makeObjFunc!Complex(64);
            Complex!F[] inps = fft.inputs!F;
            foreach(i; 0 .. 64)
                inps[i] = Complex!F(i, i);

            auto phobosResult = std.numeric.fft(inps);

            fft.fft!F();
            auto libResult = fft.outputs!F;
            import std.stdio, std.algorithm;
            foreach(i; 0 .. 64){
                assert(isClose(phobosResult[i].re, libResult[i].re, 1e-4));
                assert(isClose(phobosResult[i].im, libResult[i].im, 1e-4));
            }

            auto fft2 = makeObjFunc!std_complex_t(64);
            std_complex_t!F[] inps2 = fft2.inputs!F;
            foreach(x; 0 .. 64)
                inps2[x] = x + x*1i;

            fft2.fft!F();
            auto libResult2 = fft2.outputs!F;
            foreach(i; 0 .. 64){
                assert(isClose(phobosResult[i].re, libResult2[i].re, 1e-4));
                assert(isClose(phobosResult[i].im, libResult2[i].im, 1e-4));
            }

            auto fft3 = makeObjFunc!complex_t(64);
            complex_t!F[] inps3 = fft3.inputs!F;
            foreach(x; 0 .. 64)
                inps3[x] = x + x*1i;

            fft3.fft!F();
            auto libResult3 = fft2.outputs!F;
            foreach(i; 0 .. 64){
                assert(isClose(phobosResult[i].re, libResult3[i].re, 1e-4));
                assert(isClose(phobosResult[i].im, libResult3[i].im, 1e-4));
            }
        }
}

unittest
{
    import std.stdio;

    cfloat[] buf = new cfloat[64];
    foreach(x; 0 .. 64){
        buf[x] = x + x*1i;
        assert(buf[x] == x + x*1i);
    }

    foreach(makeObjFunc; AliasSeq!(makePhobosFFTObject, makeFFTWObject))
        foreach(F; AliasSeq!(float, double, real))
        {
            auto fft = makeObjFunc!Complex(64);
            Complex!F[] inps = fft.inputs!F;
            foreach(i; 0 .. 64)
                inps[i] = Complex!F(i, i);

            auto phobosResult = std.numeric.inverseFft(inps);

            fft.ifft!F();
            auto libResult = fft.outputs!F;
            foreach(i; 0 .. 64){
                assert(isClose(phobosResult[i].re, libResult[i].re, 1e-4));
                assert(isClose(phobosResult[i].im, libResult[i].im, 1e-4));
            }

            auto fft2 = makeObjFunc!std_complex_t(64);
            std_complex_t!F[] inps2 = fft2.inputs!F;
            foreach(x; 0 .. 64)
                inps2[x] = x + x*1i;

            fft2.ifft!F();
            auto libResult2 = fft2.outputs!F;
            foreach(i; 0 .. 64){
                assert(isClose(phobosResult[i].re, libResult2[i].re, 1e-4));
                assert(isClose(phobosResult[i].im, libResult2[i].im, 1e-4));
            }

            auto fft3 = makeObjFunc!complex_t(64);
            complex_t!F[] inps3 = fft3.inputs!F;
            foreach(x; 0 .. 64)
                inps3[x] = cpx!(complex_t, float)(x, x);

            fft3.ifft!F();
            auto libResult3 = fft3.outputs!F;
            foreach(i; 0 .. 64){
                assert(isClose(phobosResult[i].re, libResult3[i].re, 1e-4));
                assert(isClose(phobosResult[i].im, libResult3[i].im, 1e-4));
            }
        }
}

unittest
{
    import std.algorithm;

    auto fftobj = makeFFTWObject!complex_t(64);
    auto ips = fftobj.inputs!float;
    foreach(x; 0 .. 64)
        ips[x] = cpx!(complex_t, float)(x, x);

    fftobj.fft!float();

    auto fftobj2 = makeFFTWObject!complex_t(64);
    fftobj2.fftFrom!float(iota(64).map!(a => complex_t!float(a, a)));

    foreach(i; 0 .. 64){
        assert(isClose(fftobj2.outputs!float[i].re, fftobj.outputs!float[i].re));
        assert(isClose(fftobj2.outputs!float[i].im, fftobj.outputs!float[i].im));
    }
}


struct FFTObjectBank(alias makeFFTObject)
{
    alias FFTObject = typeof(makeFFTObject(0));


    ref FFTObject opIndex(size_t n)
    {
        if(n !in _bank)
            _bank[n] = makeFFTObject(n);

        return _bank[n];
    }


  private:
    static FFTObject[size_t] _bank;
}


/**

*/
deprecated
auto globalBankOf(alias makeFFTObj)() @property
{
    return FFTObjectBank!makeFFTObj.init;
}


/**
FFTを用いた離散ハートレー変換
*/
template hartleyTransform(F)
{
    void hartleyTransform(FFTObj)(ref FFTObj fftobj)
    {
        fftobj.fft!F();
        auto xs = fftobj.outputs!F;
        immutable N = xs.length;
        alias C = typeof(xs[0]);

        foreach(i; 1 .. N/2) {
            auto a = xs[i];
            auto b = xs[N-i];
            auto c = (a + b)/2;
            auto d = (a - b)/2;
            auto e = C(-d.im, d.re);

            xs[i] = c + e;
            xs[N-i] = c - e;
        }

        immutable sqrtN = sqrt(F(N));
        foreach(i; 0 .. N)
            xs[i] /= sqrtN;
    }
}

unittest
{
    float cas4(size_t i)
    {
        return cos(i*PI/2) + sin(i*PI/2);
    }

    import std.stdio;
    alias C = Complex!float;
    auto fftw = makeFFTWObject!Complex(4);
    fftw.inputs!float[0] = C(0, 0);
    fftw.inputs!float[1] = C(0, 1);
    fftw.inputs!float[2] = C(1, 0);
    fftw.inputs!float[3] = C(1, 1);
    fftw.hartleyTransform!float();
    assert(fftw.outputs!float[0].re.isClose((0 * cas4(0*0) + 0 * cas4(1*0) + 1 * cas4(2*0) + 1 * cas4(3*0))/2));
    assert(fftw.outputs!float[1].re.isClose((0 * cas4(0*1) + 0 * cas4(1*1) + 1 * cas4(2*1) + 1 * cas4(3*1))/2));
    assert(fftw.outputs!float[2].re.isClose((0 * cas4(0*2) + 0 * cas4(1*2) + 1 * cas4(2*2) + 1 * cas4(3*2))/2));
    assert(fftw.outputs!float[3].re.isClose((0 * cas4(0*3) + 0 * cas4(1*3) + 1 * cas4(2*3) + 1 * cas4(3*3))/2));

    assert(fftw.outputs!float[0].im.isClose((0 * cas4(0*0) + 1 * cas4(1*0) + 0 * cas4(2*0) + 1 * cas4(3*0))/2));
    assert(fftw.outputs!float[1].im.isClose((0 * cas4(0*1) + 1 * cas4(1*1) + 0 * cas4(2*1) + 1 * cas4(3*1))/2));
    assert(fftw.outputs!float[2].im.isClose((0 * cas4(0*2) + 1 * cas4(1*2) + 0 * cas4(2*2) + 1 * cas4(3*2))/2));
    assert(fftw.outputs!float[3].im.isClose((0 * cas4(0*3) + 1 * cas4(1*3) + 0 * cas4(2*3) + 1 * cas4(3*3))/2));

    fftw.inputs!float[] = fftw.outputs!float[];
    fftw.hartleyTransform!float();
    assert(fftw.outputs!float[0].re.isClose(0));
    assert(fftw.outputs!float[1].re.isClose(0));
    assert(fftw.outputs!float[2].re.isClose(1));
    assert(fftw.outputs!float[3].re.isClose(1));

    assert(fftw.outputs!float[0].im.isClose(0));
    assert(fftw.outputs!float[1].im.isClose(1));
    assert(fftw.outputs!float[2].im.isClose(0));
    assert(fftw.outputs!float[3].im.isClose(1));
}


/**
FFTを用いた離散ヒルベルト変換により解析信号を得ます．
fftobj.inputs!F に与えられた実数の入力信号は出力 fftobj.outputs!F の実部に保存され，出力の虚部には入力信号をヒルベルト変換した値が格納されます．
See also: https://jp.mathworks.com/help/signal/ref/hilbert.html
See also: https://ieeexplore.ieee.org/document/782222
*/
template hilbertTransform(F)
{
    void hilbertTransform(FFTObject)(ref FFTObject fftobj)
    {
        alias C = typeof(fftobj.inputs!F[0]);
        auto xs = fftobj.inputs!F;
        auto ys = fftobj.outputs!F;
        immutable N = xs.length;

        fftobj.fft!F();
        xs[0] = ys[0];
        foreach(k; 1 .. N/2) {
            xs[k] = 2 * ys[k];
        }
        xs[N/2] = ys[N/2];
        xs[N/2 + 1 .. $] = C(0);

        fftobj.ifft!F();
    }
}

unittest
{
    import std.stdio;
    alias C = Complex!float;
    auto fftw = makeFFTWObject!Complex(4);
    fftw.inputs!float[0] = C(1, 0);
    fftw.inputs!float[1] = C(2, 0);
    fftw.inputs!float[2] = C(3, 0);
    fftw.inputs!float[3] = C(4, 0);
    fftw.hilbertTransform!float();

    assert(fftw.outputs!float[0].re.isClose(1));
    assert(fftw.outputs!float[1].re.isClose(2));
    assert(fftw.outputs!float[2].re.isClose(3));
    assert(fftw.outputs!float[3].re.isClose(4));

    assert(fftw.outputs!float[0].im.isClose(1));
    assert(fftw.outputs!float[1].im.isClose(-1));
    assert(fftw.outputs!float[2].im.isClose(-1));
    assert(fftw.outputs!float[3].im.isClose(1));
}


/**
配列の前半分と後半分を入れ替えます．
この操作をFFT後の配列に適用することで，等価低域系での周波数スペクトルが得られます．
*/
void swapHalf(R)(R buf)
if(isRandomAccessRange!R && hasLength!R)
{
    import std.algorithm : swap;

    foreach(i; 0 .. buf.length/2)
        swap(buf[i], buf[$-($/2)+i]);
}

///
unittest
{
    int[] arr = [1, 2, 3, 4];
    swapHalf(arr);
    assert(arr == [3, 4, 1, 2]);

    swapHalf(arr[0 .. 3]);
    assert(arr == [1, 4, 3, 2]);
}



deprecated
Complex!F[] fftWithSwap(FftObj, F)(FftObj fftObj, in Complex!F[] buf)
in{
    assert(buf.length.isPowOf2);
}
do{
    Complex!F[] spec = new Complex!F[buf.length];
    fftObj.fft(buf, spec);
    foreach(i; 0 .. spec.length/2)
        swap(spec[i], spec[$/2+i]);

    return spec;
}


deprecated
C[] fftWithSwap(FftObj, C)(FftObj fftObj, in C[] buf)
if(is(typeof(buf[0].re) : float) && is(typeof(buf[0].im) : float))
in{
    assert(buf.length.isPowOf2);
}
do{
    alias F = typeof(buf[0].re);

    auto cpxBuf = new Complex!F[buf.length];
    foreach(i, ref e; cpxBuf)
        e = Complex!F(buf[i].re, buf[i].im);

    auto spec = fftObj.fftWithSwap(cpxBuf);
    C[] ret = new C[buf.length];

    foreach(i, ref e; ret)
        e = spec[i].re + spec[i].im * 1i;

    return ret;
}


deprecated
align(1)
struct FrequencyDomain(C)
{
    C value;
    alias value this;
}


deprecated
FrequencyDomain!C frequencyDomain(C)(C c)
{
    return typeof(return)(c);
}


deprecated
inout(C)[] convAuto(C)(inout(FrequencyDomain!C)[] array)
{
    return (cast(inout(C)*)array.ptr)[0 .. array.length];
}

//unittest
//{
//    FrequencyDomain!(Complex!float)[] carr;
//    foreach(i; 0 .. 10) carr ~= complex!float(i, 0).frequencyDomain;

//    auto carr2 = carr.convAuto!(Complex!float);
//    foreach(i; 0 .. 10){
//        assert(isClose(carr[i].re, carr2[i].re));
//        assert(isClose(carr[i].im, carr2[i].im));
//    }
//}


deprecated
cfloat[] rawReadComplex(File file, cfloat[] buf)
{
    return file.rawRead(buf);
}


deprecated
Complex!float[] rawReadComplex(File file, cfloat[] buf, Complex!float[] output)
{
    auto res = file.rawReadComplex(buf);
    foreach(i, e; res)
        output[i] = complex!float(e.re, e.im);

    return output[0 .. res.length];
}


/**
周波数領域のチャネルの推定値freqRespから時間領域のインパルス応答を推定します．
f
*/
C[] estimateImpulseResponseFromFrequencyResponse(C, Freqs)(size_t dftsize, size_t impsize, in C[] freqResp, Freqs freqIndexs, real diagonalLoading)
{
    alias R = typeof(C.init.re);
    import dffdd.utils.linalg;

    C[][] matW = new C[][](freqResp.length, impsize);
    foreach(i, f; freqIndexs.enumerate)
        foreach(j; 0 .. impsize)
            matW[i][j] = std.complex.expi(-2*PI/dftsize * f * j);

    return leastSquareEstimateTikRegSimple(matW, freqResp, diagonalLoading);
}

unittest
{
    alias C = Complex!double;

    // インパルス応答の真値
    Complex!double[] impResp = [C(1, 1), C(2, -2), C(0, 3), C(1, 0)];

    // インパルス応答から周波数応答を生成
    auto fftw = makeFFTWObject!Complex(64);
    fftw.inputs!double[0 .. 4] = impResp[0 .. 4];
    fftw.inputs!double[4 .. $] = C(0);
    fftw.fft!double();

    // 生成された周波数応答のうち，幾つかの周波数のみを抽出
    auto sampledFreqResp = fftw.outputs!double[$/8 .. $/4].dup;

    // 抽出した周波数の番号リスト
    auto freqIndexs = iota(64/8, 64/4);

    // 抽出した周波数応答からもとのインパルス応答を推定
    auto result = estimateImpulseResponseFromFrequencyResponse(64, 4, sampledFreqResp, freqIndexs, 0);
    foreach(i; 0 .. 4) {
        auto e = result[i] - impResp[i];
        assert(e.re.isClose(0, 0, 1e-6));
        assert(e.im.isClose(0, 0, 1e-6));
    }
}