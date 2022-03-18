module dffdd_simd.primitives;

public import core.simd;

version(LDC)
{
    public import ldc.attributes;

    pragma(LDC_intrinsic, "llvm.minnum.v4f32")
    float4 _mm_minps_(float4, float4) pure @safe;
    pragma(LDC_intrinsic, "llvm.maxnum.v4f32")
    float4 _mm_maxps_(float4, float4) pure @safe;
    pragma(LDC_intrinsic, "llvm.minnum.v8f32")
    float8 _mm_minps_(float8, float8) pure @safe;
    pragma(LDC_intrinsic, "llvm.maxnum.v8f32")
    float8 _mm_maxps_(float8, float8) pure @safe;
}

version(DigitalMars)
{
    float4 _mm_minps_(float4 x, float4 y) pure @safe { return cast(float4)__simd(XMM.MINPS, x, y); }
    float4 _mm_maxps_(float4 x, float4 y) pure @safe { return cast(float4)__simd(XMM.MAXPS, x, y); }
    // float8 _mm_minps_(float8 x, float8 y) pure @safe { return simd!(XMM.MINPS)(x, y); }
    // float8 _mm_maxps_(float8 x, float8 y) pure @safe { return simd!(XMM.MAXPS)(x, y); }
    enum int fastmath = 0;
}


@fastmath
V vecminmax(V, F)(V v, F vmin_, F vmax_) @trusted
{
    static if(is(V : F))
    {
        import std.algorithm : min, max;
        return max(min(v, vmax_), vmin_);
    }
    else static if(is(typeof(_mm_minps_(v, v))))
    {
        V vmin = vmin_,
          vmax = vmax_;

        return _mm_maxps_(_mm_minps_(v, vmax), vmin);
    }
    else
    {
        import std.algorithm : min, max;
        foreach(i; 0 .. V.length) {
            v[i] = max(min(v[i], vmax_), vmin_);
        }

        return v;
    }
}