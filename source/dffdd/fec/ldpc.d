module dffdd.fec.ldpc;

import std.algorithm : map, min, max;
import std.array : array;
import std.range : iota;
import std.traits : isFloatingPoint;

import dffdd.fec.fec;
import dffdd.utils.binary : Bit;
import ldpc_sp_decoder;



interface ILdpc(F) : BlockFec!F
{
    const(ubyte[][]) parityCheckMatrix() const;
    const(uint[][]) colIndexAtRow() const;
    const(uint[][]) rowIndexAtCol() const;
    uint maxRowWeight() const;
    uint maxColWeight() const;
}


class LdpcSpDecoder
{
    this(F)(in ILdpc!F ldpc, uint maxIter)
    {
        _N = ldpc.outputLength;
        _K = ldpc.inputLength;
        _row_mat = ldpc.colIndexAtRow;
        _maxIter = maxIter;

        _sp_ws_one = typeof(_sp_ws_one)(ldpc.outputLength, _row_mat);
        _sp_ws_sse = typeof(_sp_ws_sse)(ldpc.outputLength, _row_mat);
        _sp_ws_avx = typeof(_sp_ws_avx)(ldpc.outputLength, _row_mat);
    }


    Bit[] decode(F)(in F[] inputP0P1, return ref Bit[] decoded)
    {
        decoded.length = inputP0P1.length / _N * _K;

        immutable numOfBlock = decoded.length / _K;
        size_t remainBlock = numOfBlock;

        const(F)[] remainP0p1 = inputP0P1;
        Bit[] remainDecoded = decoded;

        // DMD does not support AVX version well.
        version(LDC) {
        while(remainBlock >= 8) {
            sumProductDecodeP0P1SIMD!(float[8])(_sp_ws_avx, _row_mat,  remainP0p1[0 .. _N*8], _maxIter);

            foreach(i; 0 .. 8)
                foreach(j; 0 .. _K)
                    remainDecoded[i*_K + j] = _sp_ws_avx.decoded_cw[i*_N + j];

            remainP0p1 = remainP0p1[_N*8 .. $];
            remainDecoded = remainDecoded[_K*8 .. $];
            remainBlock -= 8;
        }
        }

        while(remainBlock >= 4) {
            sumProductDecodeP0P1SIMD!(float[4])(_sp_ws_sse, _row_mat,  remainP0p1[0 .. _N*4], _maxIter);

            foreach(i; 0 .. 4)
                foreach(j; 0 .. _K)
                    remainDecoded[i*_K + j] = _sp_ws_sse.decoded_cw[i*_N + j];

            remainP0p1 = remainP0p1[_N*4 .. $];
            remainDecoded = remainDecoded[_K*4 .. $];
            remainBlock -= 4;
        }

        while(remainBlock >= 1) {
            sumProductDecodeP0P1SIMD!float(_sp_ws_one, _row_mat,  remainP0p1[0 .. _N], _maxIter);

            foreach(i; 0 .. _K)
                    remainDecoded[i] = _sp_ws_one.decoded_cw[i];

            remainP0p1 = remainP0p1[_N .. $];
            remainDecoded = remainDecoded[_K .. $];
            remainBlock -= 1;
        }

        return decoded;
    }


  private:
    size_t _N;
    size_t _K;
    const(uint[][]) _row_mat;
    uint _maxIter;
    SpDecoderWorkspace!float _sp_ws_one;
    SpDecoderWorkspace!(float[4]) _sp_ws_sse;
    SpDecoderWorkspace!(float[8]) _sp_ws_avx;
}


class BldpcIeee80211(F) : ILdpc!F
if(isFloatingPoint!F)
{
    this(uint N, uint K, uint maxIter = 20)
    {
        auto spec = getBLDPCIEEE80211Specs(N, K);
        _N = spec.N;
        _K = spec.K;
        _M = _N - _K;
        _Z = spec.Z;
        _maxIter = maxIter;

        ubyte[][] matH = new ubyte[][](_M, _N);
        foreach(ibr, brow; spec.baseH) {
            foreach(ibc, shiftval; brow) if(shiftval >= 0) {
                foreach(i; 0 .. _Z) {
                    matH[_Z*ibr + i][_Z*ibc + (i + shiftval) % _Z] = 1;
                }
            }
        }
        _matH = matH.map!"cast(immutable)a.dup".array();

        auto colIdxsOfRow = new uint[][](_M);
        auto rowIdxsOfCol = new uint[][](_N);

        foreach(ir; 0 .. _M) foreach(ic; 0 .. _N) {
            if(_matH[ir][ic] == 1) {
                rowIdxsOfCol[ic] ~= ir;
                colIdxsOfRow[ir] ~= ic;
            }
        }

        uint maxRowWeight = 0;
        foreach(ir; 0 .. _M)
            maxRowWeight = cast(uint)max(maxRowWeight, colIdxsOfRow.length);

        uint maxColWeight = 0;
        foreach(ic; 0 .. _N)
            maxColWeight = cast(uint)max(maxColWeight, rowIdxsOfCol.length);

        _colIndexAtRow = colIdxsOfRow.map!"cast(immutable)a.dup".array;
        _rowIndexAtCol = rowIdxsOfCol.map!"cast(immutable)a.dup".array;
        _maxRowWeight = maxRowWeight;
        _maxColWeight = maxColWeight;

        _decoder = new LdpcSpDecoder(this, _maxIter);
    }


    override
    Bit[] encode(in Bit[] info, return ref Bit[] encoded) const
    {
        encoded.length = info.length / _K * _N;
        foreach(i; 0 .. info.length / _K)
            encodeImpl(info[i * _K .. (i+1) * _K], encoded[i * _N .. (i+1) * _N]);

        return encoded;
    }


    override
    Bit[] encode(in Bit[] info, return ref Bit[] encoded)
    {
        return (cast(const)this).encode(info, encoded);
    }



    private
    void encodeImpl(in Bit[] info, Bit[] encoded) const
    in {
        assert(info.length == _K);
        assert(encoded.length == _N);
    }
    body {
        static struct Storage {
            size_t M;
            ubyte[] parity;
        }

        static Storage s;

        if(s.M < _M) {
            s.M = _M;
            s.parity = new ubyte[_M];
        }

        s.parity[] = 0;
        encoded[0 .. _K] = info[0 .. _K];

        // H_I u^T
        foreach(ir, colidxs; _colIndexAtRow) {
            foreach(ic; colidxs) if(ic < _K)
                s.parity[ir] += encoded[ic].value;

            s.parity[ir] %= 2;
        }

        // p_1 = ET^{-1} As + Cs
        foreach(ic; 0 .. _Z) {
            foreach(ir; iota(ic, _M, _Z))
                encoded[_K + ic].value += s.parity[ir];

            encoded[_K + ic].value %= 2;
        }

        // compute As + Bp_1
        foreach(ir, colidxs; _colIndexAtRow) {
            foreach(ic; colidxs) if(ic >= _K && ic < _K + _Z)
                s.parity[ir] += encoded[ic].value;

            s.parity[ir] %= 2;
        }

        // back substitution
        foreach(ic; iota(_K + _Z, _N, _Z)) {
            foreach(ir; 0 .. _Z) {
                encoded[ic + ir].value = s.parity[ic + ir - _K - _Z];
                s.parity[ic + ir - _K] += s.parity[ic + ir - _K - _Z];
            }
        }

        foreach(ref e; encoded)
            e.value %= 2;
    }


    override
    Bit[] decode(in F[] inputP0P1, return ref Bit[] decoded)
    {
        return _decoder.decode(inputP0P1, decoded);
    }


    bool checkParity(Bit[] codeword)
    {
        return checkCodeword(_colIndexAtRow, codeword);
    }


    override
    size_t inputLength() const { return _K; }


    override
    size_t outputLength() const { return _N; }


    override
    const(ubyte[][]) parityCheckMatrix() const
    {
        return _matH;
    }


    override
    const(uint[][]) colIndexAtRow() const
    {
        return _colIndexAtRow;
    }


    override
    const(uint[][]) rowIndexAtCol() const
    {
        return _rowIndexAtCol;
    }


    override
    uint maxRowWeight() const
    {
        return _maxRowWeight;
    }


    override
    uint maxColWeight() const
    {
        return _maxColWeight;
    }


    LdpcSpDecoder makeDecoder(int maxIter = -1) const
    {
        if(maxIter == -1)
            maxIter = this._maxIter;

        return new LdpcSpDecoder(this, maxIter);
    }


  private:
    immutable(ubyte[][]) _matH;
    immutable(uint[][]) _colIndexAtRow;
    immutable(uint[][]) _rowIndexAtCol;
    immutable uint _N, _K, _M, _Z;
    immutable uint _maxRowWeight, _maxColWeight;
    LdpcSpDecoder _decoder;
    uint _maxIter;
}

unittest
{
    import std.complex;
    import std.math;
    import std.random;
    import dffdd.mod.bpsk;
    import dffdd.mod.qpsk;
    import dffdd.mod.qam;

    alias F = float;
    alias C = Complex!F;

    BldpcIeee80211!F code = new BldpcIeee80211!F(648, 324, 20);
    
    Bit[] rndbits = [
        0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1,
        0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1,
        1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
        1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
        1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0,
        1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1,
        0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0,
        1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0,
        1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0,
        1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0,
        1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1,
        1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 0
    ].map!(a => Bit(a)).array;

    Bit[] codeword;
    code.encode(rndbits, codeword);
    assert(codeword.length == 648);

    assert(rndbits[0 .. 324] == codeword[0 .. 324]);
    assert(code.checkParity(codeword));
    assert(codeword == [
        0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1,
        0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1,
        1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
        1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
        1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0,
        1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1,
        0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0,
        1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0,
        1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 
        1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0,
        1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1,
        1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0,
        1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1,
        1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0,
        1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1,
        1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0,
        0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,
        1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1,
        0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1,
        1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1,
        0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0,
        0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0,
        1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1,
        0, 1, 0, 1, 1, 0, 1, 1
    ]);

    auto mod = BPSK!C();
    C[] syms;
    mod.modulate(codeword, syms);

    Random rnd;
    rnd.seed(0);

    // immutable N0 = 1.0/6;
    immutable N0 = 1.0;

    foreach(ref e; syms) {
        auto x = uniform01(rnd),
             y = uniform01(rnd);

        e += sqrt(-log(x)) * std.complex.expi(2*PI*y) * sqrt(N0);
    }

    auto p0p1calc = p0p1Calculator!F(mod);

    F[] p0p1;
    p0p1calc.computeP0P1(syms, p0p1, N0);

    Bit[] decoded;
    code.decode(p0p1, decoded);
    assert(decoded == rndbits);
}


struct BLDPCIEEE80211Spec
{
    immutable(int[])[] baseH;
    uint N;
    uint K;
    uint Z;
}



private
BLDPCIEEE80211Spec[size_t][size_t] constructBLDPCIEEE80211Specs()
{
    return [
    648:
        [
            /* N = 648, K = 324, R = 1/2 */
            324: BLDPCIEEE80211Spec(
                [
                    [ 0, -1, -1, -1, 0, 0, -1, -1, 0, -1, -1, 0, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 22, 0, -1, -1, 17, -1, 0, 0, 12, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 6, -1, 0, -1, 10, -1, -1, -1, 24, -1, 0, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 2, -1, -1, 0, 20, -1, -1, -1, 25, 0, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1],
                    [ 23, -1, -1, -1, 3, -1, -1, -1, 0, -1, 9, 11, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1],
                    [ 24, -1, 23, 1, 17, -1, 3, -1, 10, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 25, -1, -1, -1, 8, -1, -1, -1, 7, 18, -1, -1, 0, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 13, 24, -1, -1, 0, -1, 8, -1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 7, 20, -1, 16, 22, 10, -1, -1, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1],
                    [ 11, -1, -1, -1, 19, -1, -1, -1, 13, -1, 3, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ 25, -1, 8, -1, 23, 18, -1, 14, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 3, -1, -1, -1, 16, -1, -1, 2, 25, 5, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
                ],
                648, 324, 27
            ),

            /* N = 648, K = 432, R = 2/3 */
            432: BLDPCIEEE80211Spec(
                [
                    [ 25, 26, 14, -1, 20, -1, 2, -1, 4, -1, -1, 8, -1, 16, -1, 18, 1, 0, -1, -1, -1, -1, -1, -1],
                    [ 10, 9, 15, 11, -1, 0, -1, 1, -1, -1, 18, -1, 8, -1, 10, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 16, 2, 20, 26, 21, -1, 6, -1, 1, 26, -1, 7, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 10, 13, 5, 0, -1, 3, -1, 7, -1, -1, 26, -1, -1, 13, -1, 16, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 23, 14, 24, -1, 12, -1, 19, -1, 17, -1, -1, -1, 20, -1, 21, -1, 0, -1, -1, -1, 0, 0, -1, -1],
                    [ 6, 22, 9, 20, -1, 25, -1, 17, -1, 8, -1, 14, -1, 18, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ 14, 23, 21, 11, 20, -1, 24, -1, 18, -1, 19, -1, -1, -1, -1, 22, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 17, 11, 11, 20, -1, 21, -1, 26, -1, 3, -1, -1, 18, -1, 26, -1, 1, -1, -1, -1, -1, -1, -1, 0],
                ],
                648, 432, 27
            ),

            /* N = 648, K = 432, R = 3/4 */
            486: BLDPCIEEE80211Spec(
                [
                    [ 16, 17, 22, 24, 9, 3, 14, -1, 4, 2, 7, -1, 26, -1, 2, -1, 21, -1, 1, 0, -1, -1, -1, -1],
                    [ 25, 12, 12, 3, 3, 26, 6, 21, -1, 15, 22, -1, 15, -1, 4, -1, -1, 16, -1, 0, 0, -1, -1, -1],
                    [ 25, 18, 26, 16, 22, 23, 9, -1, 0, -1, 4, -1, 4, -1, 8, 23, 11, -1, -1, -1, 0, 0, -1, -1],
                    [ 9, 7, 0, 1, 17, -1, -1, 7, 3, -1, 3, 23, -1, 16, -1, -1, 21, -1, 0, -1, -1, 0, 0, -1],
                    [ 24, 5, 26, 7, 1, -1, -1, 15, 24, 15, -1, 8, -1, 13, -1, 13, -1, 11, -1, -1, -1, -1, 0, 0],
                    [ 2, 2, 19, 14, 24, 1, 15, 19, -1, 21, -1, 2, -1, 24, -1, 3, -1, 2, 1, -1, -1, -1, -1, 0],
                ],
                648, 432, 27
            ),

            /* N = 648, K = 540, R = 5/6 */
            540: BLDPCIEEE80211Spec(
                [
                    [ 17, 13, 8, 21, 9, 3, 18, 12, 10, 0, 4, 15, 19, 2, 5, 10, 26, 19, 13, 13, 1, 0, -1, -1],
                    [ 3, 12, 11, 14, 11, 25, 5, 18, 0, 9, 2, 26, 26, 10, 24, 7, 14, 20, 4, 2, -1, 0, 0, -1],
                    [ 22, 16, 4, 3, 10, 21, 12, 5, 21, 14, 19, 5, -1, 8, 5, 18, 11, 5, 5, 15, 0, -1, 0, 0],
                    [ 7, 7, 14, 14, 4, 16, 16, 24, 24, 10, 1, 7, 15, 6, 10, 26, 8, 18, 21, 14, 1, -1, -1, 0],
                ],
                648, 540, 27
            )
        ],

    1296:
        [
            /* N = 1296, K = 648, R = 1/2 */
            648: BLDPCIEEE80211Spec(
                [
                    [ 40, -1, -1, -1, 22, -1, 49, 23, 43, -1, -1, -1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 50, 1, -1, -1, 48, 35, -1, -1, 13, -1, 30, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 39, 50, -1, -1, 4, -1, 2, -1, -1, -1, -1, 49, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 33, -1, -1, 38, 37, -1, -1, 4, 1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1],
                    [ 45, -1, -1, -1, 0, 22, -1, -1, 20, 42, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1],
                    [ 51, -1, -1, 48, 35, -1, -1, -1, 44, -1, 18, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 47, 11, -1, -1, -1, 17, -1, -1, 51, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 5, -1, 25, -1, 6, -1, 45, -1, 13, 40, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 33, -1, -1, 34, 24, -1, -1, -1, 23, -1, -1, 46, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1],
                    [ 1, -1, 27, -1, 1, -1, -1, -1, 38, -1, 44, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ -1, 18, -1, -1, 23, -1, -1, 8, 0, 35, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 49, -1, 17, -1, 30, -1, -1, -1, 34, -1, -1, 19, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
                ],
                1296, 648, 54,
            ),

            864: BLDPCIEEE80211Spec(
                [
                    [ 39, 31, 22, 43, -1, 40, 4, -1, 11, -1, -1, 50, -1, -1, -1, 6, 1, 0, -1, -1, -1, -1, -1, -1],
                    [ 25, 52, 41, 2, 6, -1, 14, -1, 34, -1, -1, -1, 24, -1, 37, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 43, 31, 29, 0, 21, -1, 28, -1, -1, 2, -1, -1, 7, -1, 17, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 20, 33, 48, -1, 4, 13, -1, 26, -1, -1, 22, -1, -1, 46, 42, -1, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 45, 7, 18, 51, 12, 25, -1, -1, -1, 50, -1, -1, 5, -1, -1, -1, 0, -1, -1, -1, 0, 0, -1, -1],
                    [ 35, 40, 32, 16, 5, -1, -1, 18, -1, -1, 43, 51, -1, 32, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ 9, 24, 13, 22, 28, -1, -1, 37, -1, -1, 25, -1, -1, 52, -1, 13, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 32, 22, 4, 21, 16, -1, -1, -1, 27, 28, -1, 38, -1, -1, -1, 8, 1, -1, -1, -1, -1, -1, -1, 0],
                ],
                1296, 864, 54
            ),

            972: BLDPCIEEE80211Spec(
                [
                    [ 39, 40, 51, 41, 3, 29, 8, 36, -1, 14, -1, 6, -1, 33, -1, 11, -1, 4, 1, 0, -1, -1, -1, -1],
                    [ 48, 21, 47, 9, 48, 35, 51, -1, 38, -1, 28, -1, 34, -1, 50, -1, 50, -1, -1, 0, 0, -1, -1, -1],
                    [ 30, 39, 28, 42, 50, 39, 5, 17, -1, 6, -1, 18, -1, 20, -1, 15, -1, 40, -1, -1, 0, 0, -1, -1],
                    [ 29, 0, 1, 43, 36, 30, 47, -1, 49, -1, 47, -1, 3, -1, 35, -1, 34, -1, 0, -1, -1, 0, 0, -1],
                    [ 1, 32, 11, 23, 10, 44, 12, 7, -1, 48, -1, 4, -1, 9, -1, 17, -1, 16, -1, -1, -1, -1, 0, 0],
                    [ 13, 7, 15, 47, 23, 16, 47, -1, 43, -1, 29, -1, 52, -1, 2, -1, 53, -1, 1, -1, -1, -1, -1, 0],
                ],
                1296, 972, 54
            ),

            1080: BLDPCIEEE80211Spec(
                [
                    [ 48, 29, 37, 52, 2, 16, 6, 14, 53, 31, 34, 5, 18, 42, 53, 31, 45, -1, 46, 52, 1, 0, -1, -1],
                    [ 17, 4, 30, 7, 43, 11, 24, 6, 14, 21, 6, 39, 17, 40, 47, 7, 15, 41, 19, -1, -1, 0, 0, -1],
                    [ 7, 2, 51, 31, 46, 23, 16, 11, 53, 40, 10, 7, 46, 53, 33, 35, -1, 25, 35, 38, 0, -1, 0, 0],
                    [ 19, 48, 41, 1, 10, 7, 36, 47, 5, 29, 52, 52, 31, 10, 26, 6, 3, 2, -1, 51, 1, -1, -1, 0],
                ],
                1296, 1080, 54
            )
        ],

    1944:
        [
            972: BLDPCIEEE80211Spec(
                [
                    [ 57, -1, -1, -1, 50, -1, 11, -1, 50, -1, 79, -1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 3, -1, 28, -1, 0, -1, -1, -1, 55, 7, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 30, -1, -1, -1, 24, 37, -1, -1, 56, 14, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 62, 53, -1, -1, 53, -1, -1, 3, 35, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1],
                    [ 40, -1, -1, 20, 66, -1, -1, 22, 28, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1],
                    [ 0, -1, -1, -1, 8, -1, 42, -1, 50, -1, -1, 8, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 69, 79, 79, -1, -1, -1, 56, -1, 52, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 65, -1, -1, -1, 38, 57, -1, -1, 72, -1, 27, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 64, -1, -1, -1, 14, 52, -1, -1, 30, -1, -1, 32, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1],
                    [ -1, 45, -1, 70, 0, -1, -1, -1, 77, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ 2, 56, -1, 57, 35, -1, -1, -1, -1, -1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 24, -1, 61, -1, 60, -1, -1, 27, 51, -1, -1, 16, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
                ],
                1944, 972, 81
            ),

            1296: BLDPCIEEE80211Spec(
                [
                    [ 61, 75, 4, 63, 56, -1, -1, -1, -1, -1, -1, 8, -1, 2, 17, 25, 1, 0, -1, -1, -1, -1, -1, -1],
                    [ 56, 74, 77, 20, -1, -1, -1, 64, 24, 4, 67, -1, 7, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 28, 21, 68, 10, 7, 14, 65, -1, -1, -1, 23, -1, -1, -1, 75, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 48, 38, 43, 78, 76, -1, -1, -1, -1, 5, 36, -1, 15, 72, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 40, 2, 53, 25, -1, 52, 62, -1, 20, -1, -1, 44, -1, -1, -1, -1, 0, -1, -1, -1, 0, 0, -1, -1],
                    [ 69, 23, 64, 10, 22, -1, 21, -1, -1, -1, -1, -1, 68, 23, 29, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ 12, 0, 68, 20, 55, 61, -1, 40, -1, -1, -1, 52, -1, -1, -1, 44, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 58, 8, 34, 64, 78, -1, -1, 11, 78, 24, -1, -1, -1, -1, -1, 58, 1, -1, -1, -1, -1, -1, -1, 0],
                ],
                1944, 1296, 81
            ),

            1458: BLDPCIEEE80211Spec(
                [
                    [ 48, 29, 28, 39, 9, 61, -1, -1, -1, 63, 45, 80, -1, -1, -1, 37, 32, 22, 1, 0, -1, -1, -1, -1],
                    [ 4, 49, 42, 48, 11, 30, -1, -1, -1, 49, 17, 41, 37, 15, -1, 54, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 35, 76, 78, 51, 37, 35, 21, -1, 17, 64, -1, -1, -1, 59, 7, -1, -1, 32, -1, -1, 0, 0, -1, -1],
                    [ 9, 65, 44, 9, 54, 56, 73, 34, 42, -1, -1, -1, 35, -1, -1, -1, 46, 39, 0, -1, -1, 0, 0, -1],
                    [ 3, 62, 7, 80, 68, 26, -1, 80, 55, -1, 36, -1, 26, -1, 9, -1, 72, -1, -1, -1, -1, -1, 0, 0],
                    [ 26, 75, 33, 21, 69, 59, 3, 38, -1, -1, -1, 35, -1, 62, 36, 26, -1, -1, 1, -1, -1, -1, -1, 0],
                ],
                1944, 1458, 81
            ),

            1620: BLDPCIEEE80211Spec(
                [
                    [ 13, 48, 80, 66, 4, 74, 7, 30, 76, 52, 37, 60, -1, 49, 73, 31, 74, 73, 23, -1, 1, 0, -1, -1],
                    [ 69, 63, 74, 56, 64, 77, 57, 65, 6, 16, 51, -1, 64, -1, 68, 9, 48, 62, 54, 27, -1, 0, 0, -1],
                    [ 51, 15, 0, 80, 24, 25, 42, 54, 44, 71, 71, 9, 67, 35, -1, 58, -1, 29, -1, 53, 0, -1, 0, 0],
                    [ 16, 29, 36, 41, 44, 56, 59, 37, 50, 24, -1, 65, 4, 65, 52, -1, 4, -1, 73, 52, 1, -1, -1, 0],
                ],
                1944, 1620, 81
            )
        ]
    ];
}




private
BLDPCIEEE80211Spec getBLDPCIEEE80211Specs(size_t N, size_t K)
{
    if(N == 648) {
        if(K == 324){
            return BLDPCIEEE80211Spec(
                [
                    [ 0, -1, -1, -1, 0, 0, -1, -1, 0, -1, -1, 0, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 22, 0, -1, -1, 17, -1, 0, 0, 12, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 6, -1, 0, -1, 10, -1, -1, -1, 24, -1, 0, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 2, -1, -1, 0, 20, -1, -1, -1, 25, 0, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1],
                    [ 23, -1, -1, -1, 3, -1, -1, -1, 0, -1, 9, 11, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1],
                    [ 24, -1, 23, 1, 17, -1, 3, -1, 10, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 25, -1, -1, -1, 8, -1, -1, -1, 7, 18, -1, -1, 0, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 13, 24, -1, -1, 0, -1, 8, -1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 7, 20, -1, 16, 22, 10, -1, -1, 23, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1],
                    [ 11, -1, -1, -1, 19, -1, -1, -1, 13, -1, 3, 17, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ 25, -1, 8, -1, 23, 18, -1, 14, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 3, -1, -1, -1, 16, -1, -1, 2, 25, 5, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
                ],
                648, 324, 27
            );
        } else if(K == 432) {
            return BLDPCIEEE80211Spec(
                [
                    [ 25, 26, 14, -1, 20, -1, 2, -1, 4, -1, -1, 8, -1, 16, -1, 18, 1, 0, -1, -1, -1, -1, -1, -1],
                    [ 10, 9, 15, 11, -1, 0, -1, 1, -1, -1, 18, -1, 8, -1, 10, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 16, 2, 20, 26, 21, -1, 6, -1, 1, 26, -1, 7, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 10, 13, 5, 0, -1, 3, -1, 7, -1, -1, 26, -1, -1, 13, -1, 16, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 23, 14, 24, -1, 12, -1, 19, -1, 17, -1, -1, -1, 20, -1, 21, -1, 0, -1, -1, -1, 0, 0, -1, -1],
                    [ 6, 22, 9, 20, -1, 25, -1, 17, -1, 8, -1, 14, -1, 18, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ 14, 23, 21, 11, 20, -1, 24, -1, 18, -1, 19, -1, -1, -1, -1, 22, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 17, 11, 11, 20, -1, 21, -1, 26, -1, 3, -1, -1, 18, -1, 26, -1, 1, -1, -1, -1, -1, -1, -1, 0],
                ],
                648, 432, 27
            );
        } else if(K == 486) {
            return BLDPCIEEE80211Spec(
                [
                    [ 16, 17, 22, 24, 9, 3, 14, -1, 4, 2, 7, -1, 26, -1, 2, -1, 21, -1, 1, 0, -1, -1, -1, -1],
                    [ 25, 12, 12, 3, 3, 26, 6, 21, -1, 15, 22, -1, 15, -1, 4, -1, -1, 16, -1, 0, 0, -1, -1, -1],
                    [ 25, 18, 26, 16, 22, 23, 9, -1, 0, -1, 4, -1, 4, -1, 8, 23, 11, -1, -1, -1, 0, 0, -1, -1],
                    [ 9, 7, 0, 1, 17, -1, -1, 7, 3, -1, 3, 23, -1, 16, -1, -1, 21, -1, 0, -1, -1, 0, 0, -1],
                    [ 24, 5, 26, 7, 1, -1, -1, 15, 24, 15, -1, 8, -1, 13, -1, 13, -1, 11, -1, -1, -1, -1, 0, 0],
                    [ 2, 2, 19, 14, 24, 1, 15, 19, -1, 21, -1, 2, -1, 24, -1, 3, -1, 2, 1, -1, -1, -1, -1, 0],
                ],
                648, 432, 27
            );
        } else if(K == 540) {
            return BLDPCIEEE80211Spec(
                [
                    [ 17, 13, 8, 21, 9, 3, 18, 12, 10, 0, 4, 15, 19, 2, 5, 10, 26, 19, 13, 13, 1, 0, -1, -1],
                    [ 3, 12, 11, 14, 11, 25, 5, 18, 0, 9, 2, 26, 26, 10, 24, 7, 14, 20, 4, 2, -1, 0, 0, -1],
                    [ 22, 16, 4, 3, 10, 21, 12, 5, 21, 14, 19, 5, -1, 8, 5, 18, 11, 5, 5, 15, 0, -1, 0, 0],
                    [ 7, 7, 14, 14, 4, 16, 16, 24, 24, 10, 1, 7, 15, 6, 10, 26, 8, 18, 21, 14, 1, -1, -1, 0],
                ],
                648, 540, 27
            );
        }
    } else if(N == 1296) {
        if(K == 648) {
            return BLDPCIEEE80211Spec(
                [
                    [ 40, -1, -1, -1, 22, -1, 49, 23, 43, -1, -1, -1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 50, 1, -1, -1, 48, 35, -1, -1, 13, -1, 30, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 39, 50, -1, -1, 4, -1, 2, -1, -1, -1, -1, 49, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 33, -1, -1, 38, 37, -1, -1, 4, 1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1],
                    [ 45, -1, -1, -1, 0, 22, -1, -1, 20, 42, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1],
                    [ 51, -1, -1, 48, 35, -1, -1, -1, 44, -1, 18, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 47, 11, -1, -1, -1, 17, -1, -1, 51, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 5, -1, 25, -1, 6, -1, 45, -1, 13, 40, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 33, -1, -1, 34, 24, -1, -1, -1, 23, -1, -1, 46, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1],
                    [ 1, -1, 27, -1, 1, -1, -1, -1, 38, -1, 44, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ -1, 18, -1, -1, 23, -1, -1, 8, 0, 35, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 49, -1, 17, -1, 30, -1, -1, -1, 34, -1, -1, 19, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
                ],
                1296, 648, 54,
            );
        } else if(K == 864) {
            return BLDPCIEEE80211Spec(
                [
                    [ 39, 31, 22, 43, -1, 40, 4, -1, 11, -1, -1, 50, -1, -1, -1, 6, 1, 0, -1, -1, -1, -1, -1, -1],
                    [ 25, 52, 41, 2, 6, -1, 14, -1, 34, -1, -1, -1, 24, -1, 37, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 43, 31, 29, 0, 21, -1, 28, -1, -1, 2, -1, -1, 7, -1, 17, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 20, 33, 48, -1, 4, 13, -1, 26, -1, -1, 22, -1, -1, 46, 42, -1, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 45, 7, 18, 51, 12, 25, -1, -1, -1, 50, -1, -1, 5, -1, -1, -1, 0, -1, -1, -1, 0, 0, -1, -1],
                    [ 35, 40, 32, 16, 5, -1, -1, 18, -1, -1, 43, 51, -1, 32, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ 9, 24, 13, 22, 28, -1, -1, 37, -1, -1, 25, -1, -1, 52, -1, 13, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 32, 22, 4, 21, 16, -1, -1, -1, 27, 28, -1, 38, -1, -1, -1, 8, 1, -1, -1, -1, -1, -1, -1, 0],
                ],
                1296, 864, 54
            );
        } else if(K == 972) {
            return BLDPCIEEE80211Spec(
                [
                    [ 39, 40, 51, 41, 3, 29, 8, 36, -1, 14, -1, 6, -1, 33, -1, 11, -1, 4, 1, 0, -1, -1, -1, -1],
                    [ 48, 21, 47, 9, 48, 35, 51, -1, 38, -1, 28, -1, 34, -1, 50, -1, 50, -1, -1, 0, 0, -1, -1, -1],
                    [ 30, 39, 28, 42, 50, 39, 5, 17, -1, 6, -1, 18, -1, 20, -1, 15, -1, 40, -1, -1, 0, 0, -1, -1],
                    [ 29, 0, 1, 43, 36, 30, 47, -1, 49, -1, 47, -1, 3, -1, 35, -1, 34, -1, 0, -1, -1, 0, 0, -1],
                    [ 1, 32, 11, 23, 10, 44, 12, 7, -1, 48, -1, 4, -1, 9, -1, 17, -1, 16, -1, -1, -1, -1, 0, 0],
                    [ 13, 7, 15, 47, 23, 16, 47, -1, 43, -1, 29, -1, 52, -1, 2, -1, 53, -1, 1, -1, -1, -1, -1, 0],
                ],
                1296, 972, 54
            );
        } else if(K == 1080) {
            return BLDPCIEEE80211Spec(
                [
                    [ 48, 29, 37, 52, 2, 16, 6, 14, 53, 31, 34, 5, 18, 42, 53, 31, 45, -1, 46, 52, 1, 0, -1, -1],
                    [ 17, 4, 30, 7, 43, 11, 24, 6, 14, 21, 6, 39, 17, 40, 47, 7, 15, 41, 19, -1, -1, 0, 0, -1],
                    [ 7, 2, 51, 31, 46, 23, 16, 11, 53, 40, 10, 7, 46, 53, 33, 35, -1, 25, 35, 38, 0, -1, 0, 0],
                    [ 19, 48, 41, 1, 10, 7, 36, 47, 5, 29, 52, 52, 31, 10, 26, 6, 3, 2, -1, 51, 1, -1, -1, 0],
                ],
                1296, 1080, 54
            );
        }
    } else if(N == 1944) {
        if(K == 972) {
            return BLDPCIEEE80211Spec(
                [
                    [ 57, -1, -1, -1, 50, -1, 11, -1, 50, -1, 79, -1, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 3, -1, 28, -1, 0, -1, -1, -1, 55, 7, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 30, -1, -1, -1, 24, 37, -1, -1, 56, 14, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1],
                    [ 62, 53, -1, -1, 53, -1, -1, 3, 35, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1],
                    [ 40, -1, -1, 20, 66, -1, -1, 22, 28, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1],
                    [ 0, -1, -1, -1, 8, -1, 42, -1, 50, -1, -1, 8, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 69, 79, 79, -1, -1, -1, 56, -1, 52, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 65, -1, -1, -1, 38, 57, -1, -1, 72, -1, 27, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 64, -1, -1, -1, 14, 52, -1, -1, 30, -1, -1, 32, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1],
                    [ -1, 45, -1, 70, 0, -1, -1, -1, 77, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ 2, 56, -1, 57, 35, -1, -1, -1, -1, -1, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 24, -1, 61, -1, 60, -1, -1, 27, 51, -1, -1, 16, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
                ],
                1944, 972, 81
            );
        } else if(K == 1296) {
            return BLDPCIEEE80211Spec(
                [
                    [ 61, 75, 4, 63, 56, -1, -1, -1, -1, -1, -1, 8, -1, 2, 17, 25, 1, 0, -1, -1, -1, -1, -1, -1],
                    [ 56, 74, 77, 20, -1, -1, -1, 64, 24, 4, 67, -1, 7, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1],
                    [ 28, 21, 68, 10, 7, 14, 65, -1, -1, -1, 23, -1, -1, -1, 75, -1, -1, -1, 0, 0, -1, -1, -1, -1],
                    [ 48, 38, 43, 78, 76, -1, -1, -1, -1, 5, 36, -1, 15, 72, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 40, 2, 53, 25, -1, 52, 62, -1, 20, -1, -1, 44, -1, -1, -1, -1, 0, -1, -1, -1, 0, 0, -1, -1],
                    [ 69, 23, 64, 10, 22, -1, 21, -1, -1, -1, -1, -1, 68, 23, 29, -1, -1, -1, -1, -1, -1, 0, 0, -1],
                    [ 12, 0, 68, 20, 55, 61, -1, 40, -1, -1, -1, 52, -1, -1, -1, 44, -1, -1, -1, -1, -1, -1, 0, 0],
                    [ 58, 8, 34, 64, 78, -1, -1, 11, 78, 24, -1, -1, -1, -1, -1, 58, 1, -1, -1, -1, -1, -1, -1, 0],
                ],
                1944, 1296, 81
            );
        } else if(K == 1458) {
            return BLDPCIEEE80211Spec(
                [
                    [ 48, 29, 28, 39, 9, 61, -1, -1, -1, 63, 45, 80, -1, -1, -1, 37, 32, 22, 1, 0, -1, -1, -1, -1],
                    [ 4, 49, 42, 48, 11, 30, -1, -1, -1, 49, 17, 41, 37, 15, -1, 54, -1, -1, -1, 0, 0, -1, -1, -1],
                    [ 35, 76, 78, 51, 37, 35, 21, -1, 17, 64, -1, -1, -1, 59, 7, -1, -1, 32, -1, -1, 0, 0, -1, -1],
                    [ 9, 65, 44, 9, 54, 56, 73, 34, 42, -1, -1, -1, 35, -1, -1, -1, 46, 39, 0, -1, -1, 0, 0, -1],
                    [ 3, 62, 7, 80, 68, 26, -1, 80, 55, -1, 36, -1, 26, -1, 9, -1, 72, -1, -1, -1, -1, -1, 0, 0],
                    [ 26, 75, 33, 21, 69, 59, 3, 38, -1, -1, -1, 35, -1, 62, 36, 26, -1, -1, 1, -1, -1, -1, -1, 0],
                ],
                1944, 1458, 81
            );
        } else if(K == 1620) {
            return BLDPCIEEE80211Spec(
                [
                    [ 13, 48, 80, 66, 4, 74, 7, 30, 76, 52, 37, 60, -1, 49, 73, 31, 74, 73, 23, -1, 1, 0, -1, -1],
                    [ 69, 63, 74, 56, 64, 77, 57, 65, 6, 16, 51, -1, 64, -1, 68, 9, 48, 62, 54, 27, -1, 0, 0, -1],
                    [ 51, 15, 0, 80, 24, 25, 42, 54, 44, 71, 71, 9, 67, 35, -1, 58, -1, 29, -1, 53, 0, -1, 0, 0],
                    [ 16, 29, 36, 41, 44, 56, 59, 37, 50, 24, -1, 65, 4, 65, 52, -1, 4, -1, 73, 52, 1, -1, -1, 0],
                ],
                1944, 1620, 81
            );
        }
    }

    import std.exception;
    enforce(0);
    assert(0);
}
