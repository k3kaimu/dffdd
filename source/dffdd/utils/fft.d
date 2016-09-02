module dffdd.utils.fft;

import std.algorithm : swap;
import std.complex;
import std.math;
import std.numeric;
import std.stdio : File;
import std.traits : Unqual;
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
body{
    fastPODCopy(input, obj.inputs!F);
    obj.fft!F();
    fastPODCopy(obj.outputs!F, output);
}


/**

*/
void fftFrom(F, FftObj, Ri)(FftObj obj, Ri input)
if(isInputRange!Ri && isComplex!(Unqual!(ElementType!Ri), F))
in{
    static if(hasLength!Ri)
        assert(input.length == obj.inputs!F.length);
}
body{
    fastPODCopy(input, obj.inputs!F);
    obj.fft!F();
}


/**

*/
void fftTo(F, FftObj, Ro)(FftObj obj, Ro output)
if(isInputRange!Ro && hasAssignableElements!Ro)
in{
    static if(hasLength!Ro)
        assert(output.length == obj.outputs!F.length);
}
body{
    obj.fft!F();
    fastPODCopy(obj.outputs!F, output);
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
body{
    fastPODCopy(input, obj.inputs);
    obj.ifft!F();
    fastPODCopy(obj.outputs, output);
}


/**

*/
void ifftFrom(F, FftObj, Ri)(FftObj obj, Ri input)
if(isInputRange!Ri && isComplex!(Unqual!(ElementType!Ri), F))
in{
    static if(hasLength!Ri)
        assert(input.length == obj.inputs!F.length);
}
body{
    fastPODCopy(input, obj.inputs!F);
    obj.ifft!F();
}


/**

*/
void ifftTo(F, FftObj, Ro)(FftObj obj, Ro output)
if(isInputRange!Ro && hasAssignableElements!Ro)
in{
    static if(hasLength!Ro)
        assert(output.length == obj.outputs!F.length);
}
body{
    obj.ifft!F();
    fastPODCopy(obj.outputs!F, output);
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


    ~this()
    {
        if(_floatFwdPlan !is null){
            fftwf_destroy_plan(_floatFwdPlan);
            _floatFwdPlan = null;
        }

        if(_doubleFwdPlan !is null){
            fftw_destroy_plan(_doubleFwdPlan);
            _doubleFwdPlan = null;
        }

        if(_realFwdPlan !is null){
            fftwl_destroy_plan(_realFwdPlan);
            _realFwdPlan = null;
        }

        if(_floatInvPlan !is null){
            fftwf_destroy_plan(_floatInvPlan);
            _floatInvPlan = null;
        }

        if(_doubleInvPlan !is null){
            fftw_destroy_plan(_doubleInvPlan);
            _doubleInvPlan = null;
        }

        if(_realInvPlan !is null){
            fftwl_destroy_plan(_realInvPlan);
            _realInvPlan = null;
        }

        //static assert(0, "free in/out array!!!");
        if(_float_in !is null){
            fftwf_free(_float_in.ptr);
            _float_in = null;
        }

        if(_float_out !is null){
            fftwf_free(_float_out.ptr);
            _float_out = null;
        }

        if(_double_in !is null){
            fftw_free(_double_in.ptr);
            _double_in = null;
        }

        if(_double_out !is null){
            fftw_free(_double_out.ptr);
            _double_out = null;
        }

        if(_real_in !is null){
            fftwl_free(_real_in.ptr);
            _real_in = null;
        }

        if(_real_out !is null){
            fftwl_free(_real_out.ptr);
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
                fftwf_execute(_floatFwdPlan);
            else
                fftwf_execute(_floatInvPlan);
        }
        else static if(is(F == double))
        {
            static if(direction == FFTW_FORWARD)
                fftw_execute(_doubleFwdPlan);
            else
                fftw_execute(_doubleInvPlan);
        }
        else static if(is(F == real))
        {
            static if(direction == FFTW_FORWARD)
                fftwl_execute(_realFwdPlan);
            else
                fftwl_execute(_realInvPlan);
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
    auto makePlanAndArray(F)(size_t n, Cpx!F* input, Cpx!F* output)
    {
        auto newInput = input;
        auto newOutput = output;

        synchronized(fftwMutex)
        {
          static if(is(F == float))
          {
            if(newInput is null) newInput = cast(Cpx!F*)fftwf_malloc(F.sizeof * 2 * n);
            if(newOutput is null) newOutput = cast(Cpx!F*)fftwf_malloc(F.sizeof * 2 * n);
            auto planFwd = fftwf_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_FORWARD, FFTW_MEASURE);
            auto planInv = fftwf_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_BACKWARD, FFTW_MEASURE);
          }
          else static if(is(F == double))
          {
            if(newInput is null) newInput = cast(Cpx!F*)fftw_malloc(F.sizeof * 2 * n);
            if(newOutput is null) newOutput = cast(Cpx!F*)fftw_malloc(F.sizeof * 2 * n);
            auto planFwd = fftw_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_FORWARD, FFTW_MEASURE);
            auto planInv = fftw_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_BACKWARD, FFTW_MEASURE);
          }
          else static if(is(F == real))
          {
            if(newInput is null) newInput = cast(Cpx!F*)fftwl_malloc(F.sizeof * 2 * n);
            if(newOutput is null) newOutput = cast(Cpx!F*)fftwl_malloc(F.sizeof * 2 * n);
            auto planFwd = fftwl_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_FORWARD, FFTW_MEASURE);
            auto planInv = fftwl_plan_dft_1d(cast(int)n, cast(Complex!F*)newInput, cast(Complex!F*)newOutput, FFTW_BACKWARD, FFTW_MEASURE);
          }

            if(newInput is null || newOutput is null || planFwd is null || planInv is null)
            {
                import core.exception;
                onOutOfMemoryError();
            }

            return tuple(planFwd, planInv, newInput, newOutput);
        }
    }
}


/**
FFTWを利用したFFTObject
*/
alias FFTWObject(alias Cpx) = RefCounted!(FFTWObjectImpl!Cpx);


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
    foreach(F; AliasSeq!(float, double, real))
    {
        auto fft = makePhobosFFTObject!Complex(64);
        Complex!F[] inps = fft.inputs!F;
        foreach(i; 0 .. 64)
            inps[i] = Complex!F(i, i);

        auto phobosResult = std.numeric.fft(inps);

        fft.fft!F();
        auto libResult = fft.outputs!F;
        import std.stdio, std.algorithm;
        foreach(i; 0 .. 64){
            assert(approxEqual(phobosResult[i].re, libResult[i].re));
            assert(approxEqual(phobosResult[i].im, libResult[i].im));
        }

        auto fft2 = makePhobosFFTObject!std_complex_t(64);
        std_complex_t!F[] inps2 = fft2.inputs!F;
        foreach(x; 0 .. 64)
            inps2[x] = x + x*1i;

        fft2.fft!F();
        auto libResult2 = fft2.outputs!F;
        foreach(i; 0 .. 64){
            assert(approxEqual(phobosResult[i].re, libResult2[i].re));
            assert(approxEqual(phobosResult[i].im, libResult2[i].im));
        }

        auto fft3 = makePhobosFFTObject!complex_t(64);
        complex_t!F[] inps3 = fft3.inputs!F;
        foreach(x; 0 .. 64)
            inps3[x] = x + x*1i;

        fft3.fft!F();
        auto libResult3 = fft2.outputs!F;
        foreach(i; 0 .. 64){
            assert(approxEqual(phobosResult[i].re, libResult3[i].re));
            assert(approxEqual(phobosResult[i].im, libResult3[i].im));
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

    foreach(F; AliasSeq!(float, double, real))
    {
        auto fft = makeFFTWObject!Complex(64);
        Complex!F[] inps = fft.inputs!F;
        foreach(i; 0 .. 64)
            inps[i] = Complex!F(i, i);

        auto phobosResult = std.numeric.fft(inps);

        fft.fft!F();
        auto libResult = fft.outputs!F;
        foreach(i; 0 .. 64){
            assert(approxEqual(phobosResult[i].re, libResult[i].re));
            assert(approxEqual(phobosResult[i].im, libResult[i].im));
        }

        auto fft2 = makeFFTWObject!std_complex_t(64);
        std_complex_t!F[] inps2 = fft2.inputs!F;
        foreach(x; 0 .. 64)
            inps2[x] = x + x*1i;

        fft2.fft!F();
        auto libResult2 = fft2.outputs!F;
        foreach(i; 0 .. 64){
            assert(approxEqual(phobosResult[i].re, libResult2[i].re));
            assert(approxEqual(phobosResult[i].im, libResult2[i].im));
        }

        auto fft3 = makeFFTWObject!complex_t(64);
        complex_t!F[] inps3 = fft3.inputs!F;
        foreach(x; 0 .. 64)
            inps3[x] = cpx!float(x, x);

        fft3.fft!F();
        auto libResult3 = fft3.outputs!F;
        foreach(i; 0 .. 64){
            assert(approxEqual(phobosResult[i].re, libResult3[i].re));
            assert(approxEqual(phobosResult[i].im, libResult3[i].im));
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

    foreach(F; AliasSeq!(float, double, real))
    {
        auto fft = makeFFTWObject!Complex(64);
        Complex!F[] inps = fft.inputs!F;
        foreach(i; 0 .. 64)
            inps[i] = Complex!F(i, i);

        auto phobosResult = std.numeric.inverseFft(inps);

        fft.ifft!F();
        auto libResult = fft.outputs!F;
        foreach(i; 0 .. 64){
            assert(approxEqual(phobosResult[i].re, libResult[i].re));
            assert(approxEqual(phobosResult[i].im, libResult[i].im));
        }

        auto fft2 = makeFFTWObject!std_complex_t(64);
        std_complex_t!F[] inps2 = fft2.inputs!F;
        foreach(x; 0 .. 64)
            inps2[x] = x + x*1i;

        fft2.ifft!F();
        auto libResult2 = fft2.outputs!F;
        foreach(i; 0 .. 64){
            assert(approxEqual(phobosResult[i].re, libResult2[i].re));
            assert(approxEqual(phobosResult[i].im, libResult2[i].im));
        }

        auto fft3 = makeFFTWObject!complex_t(64);
        complex_t!F[] inps3 = fft3.inputs!F;
        foreach(x; 0 .. 64)
            inps3[x] = cpx!float(x, x);

        fft3.ifft!F();
        auto libResult3 = fft3.outputs!F;
        foreach(i; 0 .. 64){
            assert(approxEqual(phobosResult[i].re, libResult3[i].re));
            assert(approxEqual(phobosResult[i].im, libResult3[i].im));
        }
    }
}

unittest
{
    import std.algorithm;

    auto fftobj = makeFFTWObject!complex_t(64);
    auto ips = fftobj.inputs!float;
    foreach(x; 0 .. 64)
        ips[x] = cpx!float(x, x);

    fftobj.fft!float();

    auto fftobj2 = makeFFTWObject!complex_t(64);
    fftobj2.fftFrom!float(iota(64).map!(a => complex_t!float(a, a)));

    foreach(i; 0 .. 64){
        assert(approxEqual(fftobj2.outputs!float[i].re, fftobj.outputs!float[i].re));
        assert(approxEqual(fftobj2.outputs!float[i].im, fftobj.outputs!float[i].im));
    }
}


/**
配列の前半分と後半分を入れ替えます．
この操作をFFT後の配列に適用することで，等価低域系での周波数スペクトルが得られます．
*/
void swapHalf(C)(C[] buf)
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
C[] fftWithSwap(FftObj, C)(FftObj fftObj, in C[] buf)
in{
    assert(buf.length.isPowOf2);
}
body{
    C[] spec = new C[buf.length];
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
body{
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

unittest
{
    FrequencyDomain!(Complex!float)[] carr;
    foreach(i; 0 .. 10) carr ~= complex!float(i, 0).frequencyDomain;

    auto carr2 = carr.convAuto!(Complex!float);
    foreach(i; 0 .. 10){
        assert(approxEqual(carr[i].re, carr2[i].re));
        assert(approxEqual(carr[i].im, carr2[i].im));
    }
}


cfloat[] rawReadComplex(File file, cfloat[] buf)
{
    return file.rawRead(buf);
}


Complex!float[] rawReadComplex(File file, cfloat[] buf, Complex!float[] output)
{
    auto res = file.rawReadComplex(buf);
    foreach(i, e; res)
        output[i] = complex!float(e.re, e.im);

    return output[0 .. res.length];
}


