module dffdd.filter.orthogonalize;

import carbon.math;
//import carbon.linear;
import dffdd.utils.linalg;
import dffdd.filter.traits;

import std.math;
import std.range;
import std.complex;
import std.traits;
import std.meta;

import mir.ndslice : Slice, slice, Contiguous, sliced, diagonal;

final class DiagonalizationOBFFactory(C, basisFuncs...)
{
    enum size_t Dim = basisFuncs.length;
    alias F = typeof(C.init.re);

  static if(is(F == float))
    alias BuiltInCpx = cfloat;
  else
    alias BuiltInCpx = cdouble;

    this()
    {
        _covm = (new C[Dim * Dim]).sliced(Dim, Dim);
    }


    void start()
    {
        _cnt = 0;
        foreach(i; 0 .. Dim) foreach(j; 0 .. Dim)
            _covm[i, j] = complexZero!C;
    }


    void finish()
    {
        import std.stdio;
        _covm[] /= _cnt;

        F[] eigenvalues = new F[Dim];
        heev(Dim, cast(BuiltInCpx*)&(_covm[0, 0]), eigenvalues.ptr);

        foreach(i; 0 .. Dim){
            eigenvalues[i] = sqrt(eigenvalues[i]);
            foreach(j; 0 .. Dim){
                _covm[i, j] = _covm[i, j].conj;
                _covm[i, j] /= eigenvalues[i];
            }
        }
    }


    void put(R)(R xs)
    if(isInputRange!R)
    {
        foreach(x; xs){
            C[Dim] fs;
            foreach(i, f; basisFuncs){
                fs[i] = f(x);
            }

            foreach(i; 0 .. Dim)
                foreach(j; 0 .. Dim){
                    _covm[i, j] += fs[i] * fs[j].conj;
                }
            ++_cnt;
        }
    }


    void getCoefs(C1)(size_t idx, C1[] buf)
    {
        foreach(i; 0 .. Dim)
            buf[i] = _covm[idx, i];
    }


  private:
    Slice!(Contiguous, 2, C*) _covm;
    size_t _cnt;
}


final class GramSchmidtOBFFactory(C)
{
    alias F = typeof(C.init.re);

    this(size_t dim)
    {
        _m = new C[dim * dim].sliced(dim, dim);
        _dim = dim;
    }


    void start()
    {
        foreach(i; 0 .. _dim) foreach(j; 0 .. _dim)
            _m[i, j] = complexZero!C;

        _xs.length = 0;
    }


    void finish()
    {
        // 行iとjの内積
        C innerProduct(size_t i, size_t j)
        {
            C sum = complexZero!C;
            foreach(es; _xs){
                C v1 = complexZero!C,
                  v2 = complexZero!C;
                foreach(k, e; es){
                    v1 += _m[i, k] * e;
                    v2 += _m[j, k] * e;
                }

                sum += v1 * v2.conj;
            }

            sum /= _xs.length;
            return sum;
        }


        // 
        void normalize(size_t i)
        {
            F sum = 0;
            foreach(es; _xs){
                C v1 = complexZero!C;
                foreach(k, e; es)
                    v1 += _m[i, k] * e;

                sum += v1.re^^2 + v1.im^^2;
            }

            sum /= _xs.length;
            sum = sqrt(sum);
            foreach(k; 0 .. _dim)
                _m[i, k] /= sum;
        }

        import std.stdio;
        import std.conv;
        import std.exception;

        auto res = new C[_dim];
        foreach(i; 0 .. _dim){
            _m[i, i] = complexZero!C + 1;

            foreach(j; 0 .. _dim) res[j] = complexZero!C;
            foreach(j; 0 .. i){
                immutable c = innerProduct(i, j);
                std.exception.enforce(!isNaN(c.re^^2 + c.im^^2), std.conv.to!string(i) ~ std.conv.to!string(j) ~ std.conv.to!string(_m));
                foreach(k; 0 .. j+1)
                    res[k] += c * _m[j, k];
            }
            foreach(j; 0 .. i)
                _m[i, j] = -res[j];

            normalize(i);
        }
    }


    void put(R)(R xss)
    if(isInputRange!R)
    {
        foreach(xs; xss){
            assert(xs.length == _dim);
            _xs ~= xs.dup;
        }
    }


    void getCoefs(size_t idx, C[] buf)
    {
        foreach(i; 0 .. _dim)
            buf[i] = _m[idx, i];
    }


    Slice!(C*, 2, Contiguous) convertMatrix() @property
    {
        return _m;
    }


  private:
    size_t _dim;
    Slice!(C*, 2, Contiguous) _m;
    immutable(C[])[] _xs;
}


final class VectorConverter(C)
{
    alias OutputElementType(A) = C;


    this(size_t dim)
    {
        this(new C[dim * dim].sliced(dim, dim));
        this.convertMatrix[] = complexZero!C;
        this.convertMatrix.diagonal[] = 1;
    }


    this(Slice!(C*, 2, Contiguous) convMatrix)
    {
        _m = convMatrix;
        _dim = convMatrix.length!0;
    }


    Slice!(C*, 2, Contiguous) convertMatrix() @property 
    {
        return _m;
    }


    void opCallImpl(in C[] input, ref C[] output)
    {
        output.length = _dim;
        foreach(p; 0 .. _dim){
            output[p] = 0;
            foreach(q; 0 .. _dim)
                output[p] += input[q] * _m[p, q];
        }
    }

    mixin ConverterOpCalls!(const(C[]), C[]);


  private:
    size_t _dim;
    Slice!(C*, 2, Contiguous) _m;
}

unittest
{
    import std.math;

    alias C = Complex!float;

    auto conv = new VectorConverter!(Complex!float)(2);
    conv.convertMatrix[] = complexZero!C;
    conv.convertMatrix[0, 0] = 1;
    conv.convertMatrix[1, 1] = 1;

    auto res = conv([Complex!float(1, 1), Complex!float(0, 1)]);
    assert(approxEqual(res[0].re, 1));
    assert(approxEqual(res[0].im, 1));
    assert(approxEqual(res[1].re, 0));
    assert(approxEqual(res[1].im, 1));

    conv.convertMatrix[0, 1] = 1;
    res = conv([Complex!float(1, 1), Complex!float(0, 1)]);
    assert(approxEqual(res[0].re, 1));
    assert(approxEqual(res[0].im, 2));
    assert(approxEqual(res[1].re, 0));
    assert(approxEqual(res[1].im, 1));

    conv.convertMatrix[1, 0] = -1;
    res = conv([Complex!float(1, 1), Complex!float(0, 1)]);
    assert(approxEqual(res[0].re, 1));
    assert(approxEqual(res[0].im, 2));
    assert(approxEqual(res[1].re, -1));
    assert(approxEqual(res[1].im, 0));
}


final class ConvertedVectorDistorter(C, Distorter, Converter)
{
    // enum size_t dim = Distorter.dim;
    size_t outputDim() const @property { return _distorter.outputDim; }


    this(Distorter distorter, Converter converter)
    {
        _distorter = distorter;
        _converter = converter;
    }


    Distorter distorter() @property { return _distorter; }
    Converter converter() @property { return _converter; }


    void opCallImpl(C c, ref C[] output)
    {
        _converter(_distorter(c), output);
    }


    mixin ConverterOpCalls!(const(C), C[]);

  private:
    Distorter _distorter;
    Converter _converter;
}

unittest
{
    import std.math;
    import dffdd.filter.primitives;

    alias C = Complex!float;

    auto conv = new VectorConverter!(Complex!float)(2);
    auto dist = new Distorter!(Complex!float, x => x, x => x*2)();

    conv.convertMatrix[] = complexZero!C;
    conv.convertMatrix[0, 0] = 1;
    conv.convertMatrix[0, 1] = -1;
    conv.convertMatrix[1, 0] = -1;
    conv.convertMatrix[1, 1] = -1;

    auto conn = new ConvertedVectorDistorter!(Complex!float, typeof(dist), typeof(conv))(dist, conv);

    auto res = conn(Complex!float(1, 2));
    assert(approxEqual(res[0].re, -1));
    assert(approxEqual(res[0].im, -2));
    assert(approxEqual(res[1].re, -3));
    assert(approxEqual(res[1].im, -6));
}


final class OrthogonalizedVectorDistorter(C, Distorter, Orthogonalizer)
if(Distorter.inputBlockLength == 1)
{
    this(Distorter distorter, Orthogonalizer orthogonalizer)
    {
        _dist = distorter;
        _orth = orthogonalizer;
        _conn = new typeof(_conn)(_dist, new VectorConverter!C(_dist.outputDim));
        _distbuf = new C[_dist.outputDim];
    }


  static if(is(typeof((){ enum x = Distorter.outputDim; })))
    enum size_t outputDim = Distorter.outputDim;
  else
    size_t outputDim() const @property { return _dist.outputDim; }

    enum size_t inputBlockLength = 1;


    void start()
    {
        _orth.start();
    }


    void finish()
    {
        _orth.finish();
        _conn.converter.convertMatrix[] = _orth.convertMatrix;
    }


    void put(C e)
    {
        _dist(e, _distbuf);
        .put(_orth, _distbuf);
    }


    void learn(R)(R input)
    if(isInputRange!R && is(Unqual!(ElementType!R) : C))
    {
        this.start();

        foreach(e; input){
            this.put(e);
        }

        this.finish();
    }


    void opCallImpl(C input, ref C[] output)
    {
        _conn(input, output);
    }


    mixin ConverterOpCalls!(const(C), C[]);


    Distorter distorter() @property { return _dist; }
    Orthogonalizer orthogonalizer() @property { return _orth; }
    VectorConverter!C converter() @property { return _conn.converter; }
    ConvertedVectorDistorter!(C, Distorter, VectorConverter!C) convertedDistorter() @property { return _conn; }


  private:
    Distorter _dist;
    Orthogonalizer _orth;
    ConvertedVectorDistorter!(C, Distorter, VectorConverter!C) _conn;
    C[] _distbuf;
}

unittest
{
    import std.stdio;
    import dffdd.filter.primitives;
    import std.algorithm;

    alias C = Complex!float;
    alias Dist = Distorter!(Complex!float, x => x, x => x*2+1);
    alias Gs = GramSchmidtOBFFactory!C;

    auto orth = new OrthogonalizedVectorDistorter!(C, Dist, Gs)(new Dist(), new Gs(2));
    static assert(isBlockConverter!(typeof(orth), C, C[]));


    auto input = std.range.iota(4).map!(a => cast(C)std.complex.expi(PI_2 * a)).array();

    orth.learn(input);

    auto cm = orth._conn.converter.convertMatrix;

    assert(approxEqual(cm[0, 0].re, 1));
    assert(approxEqual(cm[0, 1].re, 0));
    assert(approxEqual(cm[1, 0].re, -2));
    assert(approxEqual(cm[1, 1].re, 1));

    assert(approxEqual(cm[0, 0].im, 0));
    assert(approxEqual(cm[0, 1].im, 0));
    assert(approxEqual(cm[1, 0].im, 0));
    assert(approxEqual(cm[1, 1].im, 0));

    auto output = orth(input);

    //assert(output[0].re == input)
    //
    foreach(e; zip(input, output)){
        assert(approxEqual(e[1][0].re, e[0].re));
        assert(approxEqual(e[1][0].im, e[0].im));

        assert(approxEqual(e[1][1].re, 1));
        assert(approxEqual(e[1][1].im, 0));
    }
}


final class OrthogonalizedSpectrumConverter(C, Orthogonalizer)
{
    this(Orthogonalizer orthogonalizer, size_t dim, size_t nFFT)
    {
        _orth = orthogonalizer;
        _dim = dim;

        foreach(i; 0 .. nFFT)
            _convs ~= new VectorConverter!C(dim);
    }


    size_t inputBlockLength() const @property { return _convs.length; }


    void opCallImpl(C[][] input, ref C[][] output)
    in{
        foreach(e; input) assert(e.length == _dim);
        // foreach(e; output) assert(e.length == _dim);
    }
    body{
        output.length = input.length;
        foreach(i; 0 .. _convs.length)
            _convs[i](input[i], output[i]);
    }


    mixin ConverterOpCalls!(C[][], C[][]);


    void learn(R)(R input, in bool[] scMap)
    if(isInputRange!R && is(Unqual!(ElementType!R) : const(C[][])))
    {
        const(C[][])[] specs;
        foreach(e; input)
            specs ~= e;

        foreach(f; 0 .. _convs.length){
            if(scMap[f]){
                _orth.start();
                foreach(e; specs)
                    .put(_orth, e[f]);
                _orth.finish();
                _convs[f].convertMatrix[] = _orth.convertMatrix;
            }
        }
    }


    OrthogonalizedSpectrumConverter inverter() @property
    {
        immutable nFFT = _convs.length;
        typeof(return) dst = new OrthogonalizedSpectrumConverter(_orth, _dim, nFFT);

        import std.stdio;
        writeln(nFFT);
        foreach(i; 0 .. nFFT)
            dst._convs[i].convertMatrix[] = this._convs[i].convertMatrix.transposed;
        
        return dst;
    }


    VectorConverter!C converter(size_t i) { return _convs[i]; }


  private:
    size_t _dim;
    Orthogonalizer _orth;
    VectorConverter!C[] _convs;
}


template OBFEval(basisFuncs...)
{
    C OBFEval(C)(C x, C[] coefs)
    {
        C ret = complexZero!C;
        foreach(i, bf; basisFuncs){
            ret += bf(x) * coefs[i];
        }

        return ret;
    }
}


unittest
{
    /**
    <v, w> := E[vw*],
    ||v|| := sqrt(<v, v>) = sqrt(E[xx*])
    x[n] = [-1, 0, 1]

    w0 = 1,
    u0 = w0/||w0|| = 1,
    
    w1 = x - 1<w1, u0> = x - E[x] = x,
    u1 = w1/||w1|| = w1/E[xx*] = w1/sqrt((1+0+1)/3) = sqrt(1.5)*x,
    
    w2 = x^2 - 1<w2, u0> - x<w2, u1> = x^2 - E[x^2] - x*E[x^2 * sqrt(1.5)x*]
                                     = x^2 - (2/3) - sqrt(1.5)*x*E[x * |x|^2]
                                     = x^2 - (2/3)
    u2 = w2/||w2|| = w2 / sqrt(E[(x^2 - (2/3))((x^2)* - (2/3))])
                   = w2 / sqrt(E[|x|^4 - (2/3)x^2 - (2/3)(x^2)* + 4/9])
                   = w2 / sqrt(E[|x|^4] - (2/3)E[x^2] - (2/3)E[(x^2)*] + 4/9)
                   = w2 / sqrt(2/3 - (2/3)*(2/3) - (2/3)(2/3) + 4/9)
                   = w2 / sqrt(6/9 - 4/9)
                   = w2 / sqrt(2/9)
                   = 3*w2 / sqrt(2)

    */
}
