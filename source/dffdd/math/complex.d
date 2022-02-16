module dffdd.math.complex;


import std.meta;
import std.complex;
import std.traits;

import mir.complex;


alias MirComplex = mir.complex.Complex;
alias StdComplex = std.complex.Complex;
alias mirComplex = mir.complex.complex;
alias stdComplex = std.complex.complex;

alias Complex = std.complex.Complex;
alias complex = std.complex.complex;


enum isComplex(T) = !isIntegral!T && !isFloatingPoint!T && is(typeof((T v){ auto a = v.re; auto b = v.im; }));
enum isNarrowComplex(T) = isComplex!T && is(typeof(T.init.tupleof) == std.meta.AliasSeq!(typeof(T.init.re), typeof(T.init.re)))
                        && (Unqual!T t){
                            t = T(cast(typeof(T.init.re)) 0, cast(typeof(T.init.re)) 1);
                            if(!(t.re == 0 && t.im == 1))   return false;
                            if(t.tupleof != AliasSeq!(cast(typeof(T.init.re)) 0, cast(typeof(T.init.re)) 1))
                                return false;

                            return true;
                        }(T.init);

enum isComplex(T, F) = isComplex!T && is(typeof(T.init.re) == F);


unittest
{
    static assert(isComplex!(MirComplex!float));
    static assert(isComplex!(MirComplex!double));
    static assert(isComplex!(StdComplex!float));
    static assert(isComplex!(StdComplex!double));
    static assert(!isComplex!(float));
    static assert(!isComplex!(double));

    static assert(isNarrowComplex!(MirComplex!float));
    static assert(isNarrowComplex!(MirComplex!double));
    static assert(isNarrowComplex!(StdComplex!float));
    static assert(isNarrowComplex!(StdComplex!double));
    static assert(!isNarrowComplex!(float));
    static assert(!isNarrowComplex!(double));

    static struct Cpx { float re; float im; float dummy; }
    static assert(!isNarrowComplex!Cpx);
}


inout(MirComplex!(typeof(C.init.re)))[] toMirComplexArray(C)(inout(C)[] input) @trusted
if(isNarrowComplex!C)
{
    return (cast(inout(MirComplex!(typeof(C.init.re)))*) input.ptr)[0 .. input.length];
}


inout(StdComplex!(typeof(C.init.re)))[] toStdComplexArray(C)(inout(C)[] input) @trusted
if(isNarrowComplex!C)
{
    return (cast(inout(StdComplex!(typeof(C.init.re)))) input.ptr)[0 .. input.length];
}


alias RealPartType(C) = typeof(C.init.re);



RealPartType!C abs(C)(C x)
if(isComplex!C)
{
    return std.complex.abs(StdComplex!(RealPartType!C)(x.re, x.im));
}


RealPartType!C sqAbs(C)(C x)
if(isComplex!C)
{
    return std.complex.sqAbs(StdComplex!(RealPartType!C)(x.re, x.im));
}


RealPartType!C arg(C)(C x)
if(isComplex!C)
{
    return std.complex.arg(StdComplex!(RealPartType!C)(x.re, x.im));
}

RealPartType!C norm(C)(C x)
if(isComplex!C)
{
    return std.complex.norm(StdComplex!(RealPartType!C)(x.re, x.im));
}


C conj(C)(C x)
if(isNarrowComplex!C)
{
    return C(x.re, -x.im);
}


C sin(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.sin(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C cos(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.cos(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C tan(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.tan(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C asin(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.asin(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C acos(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.acos(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C atan(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.atan(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C sinh(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.sinh(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C cosh(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.cosh(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C tanh(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.tanh(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C asinh(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.asinh(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C acosh(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.acosh(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C atanh(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.atanh(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C expi(C)(real x)
if(isNarrowComplex!C)
{
    auto y = std.complex.expi(x);
    return C(y.re, y.im);
}


C coshisinh(C)(real x)
if(isNarrowComplex!C)
{
    auto y = std.complex.coshisinh(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C sqrt(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.sqrt(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C exp(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.exp(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C log(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.log(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C log10(C)(C x)
if(isNarrowComplex!C)
{
    auto y = std.complex.log10(StdComplex!(RealPartType!C)(x.re, x.im));
    return C(y.re, y.im);
}


C pow(C1, C2)(C1 x, C2 n)
if(isNarrowComplex!C1 && (isComplex!C2 || isFloatingPoing!C2 || isInteger!C2))
{
    static if(isFloatingPoing!C2 || isInteger!C2) {
        auto y = std.complex.pow(StdComplex!(RealPartType!C)(x.re, x.im), n);
        return C(y.re, y.im);
    } else {
        auto y = std.complex.pow(StdComplex!(RealPartType!C)(x.re, x.im), StdComplex!(RealPartType!C)(n.re, n.im));
        return C(y.re, y.im);
    }
}
