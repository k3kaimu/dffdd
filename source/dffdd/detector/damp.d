module dffdd.detector.damp;

import std.math;

import dffdd.detector.primitives;

import dffdd.damp.damp;
import dffdd.math.complex;
import dffdd.math.matrix;
import dffdd.math.vector;
import dffdd.math.linalg;

import dffdd.mod.qpsk;
import dffdd.mod.qam;

import mir.ndslice : slice, sliced, Contiguous, SliceKind;


class DAMPDectector(C, Mod) : IDetector!(C, C)
if(is(Mod == QPSK!C) || is(Mod == QAM!C) && isComplex!C && (is(typeof(C.init.re) == float) || is(typeof(C.init.re) == double)) )
{
    alias F = typeof(C.init.re);

    alias InputElementType = C;
    alias OutputElementType = C;


    this(Mat)(Mod mod, in Mat chMat, size_t maxIter = 20, F damp1_1 = 1, F damp1_2 = 1, F damp2_1 = 1, F damp2_2 = 1)
    if(isMatrixLike!Mat)
    {
        _mod = mod;
        _chMat = matrix!F(chMat.length!0 * 2, chMat.length!1 * 2);
        _rowScales = vector!F(chMat.length!0 * 2);
        _chMat[] = chMat.toRealRIIR;

        // 大システム極限での正規化条件（lim |a_n|^2 = 1）
        foreach(i; 0 .. chMat.length!0 * 2) {
            auto rvec = _chMat.rowVec(i);
            auto p = 1/sqrt(norm2(rvec));
            _rowScales[i] = p;
            rvec.sliced()[] *= p;
        }

        _delta = _chMat.length!0 * 1.0 / _chMat.length!1;
        _maxIter = maxIter;
        _damp1_1 = damp1_1;
        _damp1_2 = damp1_2;
        _damp2_1 = damp2_1;
        _damp2_2 = damp2_2;
    }


    size_t inputLength() const { return _chMat.length!0 / 2; }
    size_t outputLength() const { return _chMat.length!1 / 2; }



    C[] detect(in C[] received, return ref C[] detected) const
    in(received.length % this.inputLength == 0)
    {
        immutable M = this.inputLength;
        immutable N = this.outputLength;

        if(detected.length != received.length / M * N)
            detected.length = received.length / M * N;

        // auto reals = received.toRealRI;
        auto inpvecs = vector!F(M * 2);
        foreach(n; 0 .. received.length / M) {
            inpvecs[] = received[M * n  .. M * (n+1)].sliced.vectored.toRealRI;
            foreach(i; 0 .. inpvecs.length)
                inpvecs[i] *= _rowScales[i];

            // import std.stdio;
            // writeln(inpvecs.sliced);
            Vector!(F, Contiguous) res;

            switch(_mod.symInputLength) {
                case 2: // QPSK
                    res = BODAMP!(softThrBayesOpt2PAM_QPSK!F, F)
                        (inpvecs, _chMat, _delta, 0.5, _maxIter, _damp1_1, _damp1_2, _damp2_1, _damp2_2);
                    break;
                case 4: // 16QAM
                    res = BODAMP!(softThrBayesOpt4PAM_16QAM!F, F)
                        (inpvecs, _chMat, _delta, 0.5, _maxIter, _damp1_1, _damp1_2, _damp2_1, _damp2_2);
                    break;
                case 6: // 64QAM
                    res = BODAMP!(softThrBayesOpt8PAM_64QAM!F, F)
                        (inpvecs, _chMat, _delta, 0.5, _maxIter, _damp1_1, _damp1_2, _damp2_1, _damp2_2);
                    break;
                default: {
                    import std.conv;
                    assert(0, "Unsupported modulation scheme: " ~ _mod.to!string);
                }
            }
            // writeln(res);
            detected[N * n .. N * (n+1)].sliced.vectored[] = res.fromRealRI;
        }

        return detected;
    }


//   private:
    Mod _mod; 
    Matrix!(F, Contiguous) _chMat;
    Vector!(F, Contiguous) _rowScales;
    double _delta;
    double _damp1_1;
    double _damp1_2;
    double _damp2_1;
    double _damp2_2;
    size_t _maxIter;
}

unittest
{
    import std.stdio;
    import mir.ndslice : sliced;
    import dffdd.math.exprtemplate;

    alias C = MirComplex!float;
    auto chMat = [C(1, 0), C(0, 0), C(0, 0), C(1, 0)].sliced(2, 2).matrixed;
    auto recv = [C(1/SQRT2, 1/SQRT2), C(-1/SQRT2, 1/SQRT2)];

    auto detector = new DAMPDectector!(C, QPSK!C)(QPSK!C(), chMat, 20);
    C[] dst;
    dst = detector.detect(recv, dst);
    assert(dst[0].re.isClose(recv[0].re));
    assert(dst[0].im.isClose(recv[0].im));
    assert(dst[1].re.isClose(recv[1].re));
    assert(dst[1].im.isClose(recv[1].im));
    
}


// unittest
// {
//     import dffdd.math.exprtemplate;
//     import std.stdio;
//     import mir.ndslice;

//     alias F = float;

//     auto chMat = matrix!F(4, 4);
//     chMat[] = ConstEye!(F, 1)(4);
//     writeln(chMat);

//     enum F[] arrR = [-SQRT1_2, SQRT1_2];

//     auto vec = (cast(F[])[SQRT1_2, -SQRT1_2, SQRT1_2, SQRT1_2]).sliced.vectored;
//     // pragma(msg, typeof(vec));
//     // auto res = BODAMP!(softThrBayesOpt!([0.5, 0.5], arrR, F), F)(vec, chMat, 1, 0.5, 50);
//     auto res = BODAMP!(softThrBayesOpt2PAM_QPSK!F, F)(vec, chMat, 1, 0.5, 50);

//     writeln(res);
//     writeln(makeViewOrNewSlice(vec - chMat * res, matvecAllocator));
// }
