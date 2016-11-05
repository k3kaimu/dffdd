module dffdd.filter.orthogonalize;

import carbon.math;
import carbon.linear;
import dffdd.utils.linalg;

import std.math;
import std.range;
import std.complex;
import std.experimental.ndslice;

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
    Slice!(2, C*) _covm;
    size_t _cnt;
}


final class GramSchmidtOBFFactory(C, basisFuncs...)
{
    enum size_t Dim = basisFuncs.length;
    alias F = typeof(C.init.re);

    this()
    {
        _m = new C[Dim * Dim].sliced(Dim, Dim);
    }


    void start()
    {
        foreach(i; 0 .. Dim) foreach(j; 0 .. Dim)
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
            foreach(k; 0 .. Dim)
                _m[i, k] /= sum;
        }

        import std.stdio;
        import std.conv;
        import std.exception;

        foreach(i; 0 .. Dim){
            _m[i, i] = complexZero!C + 1;

            C[Dim] res;
            foreach(j; 0 .. Dim) res[j] = complexZero!C;
            foreach(j; 0 .. i){
                auto c = innerProduct(i, j);
                std.exception.enforce(!isNaN(c.re^^2 + c.im^^2), std.conv.to!string(i) ~ std.conv.to!string(j) ~ std.conv.to!string(_m));
                foreach(k; 0 .. j+1)
                    res[k] += c * _m[j, k];
            }
            foreach(j; 0 .. i)
                _m[i, j] = -res[j];

            normalize(i);
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

            _xs ~= fs;
        }
    }


    void getCoefs(size_t idx, C[] buf)
    {
        foreach(i; 0 .. Dim)
            buf[i] = _m[idx, i];
    }


  private:
    Slice!(2, C*) _m;
    C[Dim][] _xs;
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
