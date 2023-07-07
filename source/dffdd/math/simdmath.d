module dffdd.math.simdmath;

import core.simd;
// version(D_AVX):

import inteli.avxintrin;


version(LDC)
{
    import ldc.attributes : fastmath;
}
else
{
    enum fastmath = 0;
}



float8 min(float8 x, float8 y) pure nothrow @safe @fastmath
{
    return _mm256_min_ps(x, y);
}


/**
From: https://gist.github.com/jrade/293a73f89dfef51da6522428c857802d
*/
float8 approx_exp(float8 x) nothrow @safe @fastmath
{
    enum float a = (1 << 23) / 0.69314718f;
    enum float b = (1 << 23) * (127 - 0.043677448f);
    x = a * x + b;

    enum float c = (1 << 23);
    enum float d = (1 << 23) * 255;

    immutable mask_c = _mm256_cmp_ps!_CMP_GT_OS(c, x);           // (c > x) ? 0xFFFFFFFF : 0x0
    immutable blend_c = _mm256_blendv_ps(x, 0f, mask_c);         // (c > x) ? 0 : x
    immutable mask_d = _mm256_cmp_ps!_CMP_GT_OS(x, d);           // (x > d) ? 0xFFFFFFFF : 0x0
    immutable blend_d = _mm256_blendv_ps(blend_c, d, mask_d);    // (x > d) ? d : x

    immutable n = _mm256_cvttps_epi32(blend_d);
    return _mm256_castsi256_ps(n);
}
