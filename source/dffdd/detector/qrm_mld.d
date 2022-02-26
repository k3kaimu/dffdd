module dffdd.detector.qrm_mld;

import dffdd.math.complex;
import dffdd.math.linalg;
import dffdd.math.matrix;
import dffdd.math.vector;

import mir.ndslice : Contiguous, sliced;

import dffdd.mod.primitives;
import dffdd.detector.primitives;
import std.typecons;

struct SphereDetector(C)
{
    this(const(C)[] points, F radius, Flag!"dynamicUpdate" dynamicUpdate = Yes.dynamicUpdate)
    {
        _points = points;
        _radius = radius;
        _dynamic = dynamicUpdate;
    }


    F detect(Vector!(C, Contiguous) recvY, Matrix!(C, Contiguous) chMatR, C[] decoded)
    in(decoded.length == chMatR.length!1)
    {
        _recvY = recvY;
        _chMatR = chMatR;
        _minError = F.infinity;
        _minVector = vector!C(decoded.length);

        foreach(p; _points) {
            decoded[$-1] = p;
            detectImpl(0, decoded.sliced.vectored, chMatR.length!0 - 1);
        }

        decoded.sliced.vectored()[] = _minVector;
        return _minError;
    }


    void detectImpl(F sumerror, Vector!(C, Contiguous) decoded, size_t row)
    {
        F error = sumerror + sqAbs(_recvY[row] - dot(_chMatR.rowVec(row), decoded));

        if(error >= this._radius) {
            return;
        }

        if(_dynamic && error >= _minError) {
            return;
        }

        if(row == 0) {
            if(error < _minError) {
                _minError = error;
                _minVector[] = decoded;
            }
        } else {
            foreach(p; _points) {
                decoded.sliced()[0 .. row] = C(0);
                decoded[row - 1] = p;
                detectImpl(error, decoded, row - 1);
            }
        }
    }


  private:
    alias F = typeof(C.init.re);

    const(C)[] _points;
    F _radius;
    bool _dynamic;
    Vector!(C, Contiguous) _recvY;
    Matrix!(C, Contiguous) _chMatR;
    F _minError;
    Vector!(C, Contiguous) _minVector;
}


struct MAlgorithmDetector(C)
{
    this(const(C)[] points, size_t selectM)
    {
        _points = points;
        _selectM = selectM;
    }


    F detect(Vector!(C, Contiguous) recvY, Matrix!(C, Contiguous) chMatR, C[] decoded)
    {
        import std.algorithm : sort, topNCopy, min, reduce;

        // immutable greedyError = detectImplGreedy(decoded);
        // import std.typecons;

        State[] states = [ State(0, vector!C(chMatR.length!1, C(0, 0))) ];
        State[] newstate;
        foreach_reverse(row; 0 .. chMatR.length!0)
        {
            foreach(s; states) {
                foreach(p; _points) {
                    auto newvec = vector!C(chMatR.length!1, C(0, 0));
                    newvec[] = s.vec;
                    newvec[row] = p;
                    F newerror = s.error + sqAbs(recvY[row] - dot(chMatR.rowVec(row), newvec));

                    // if(newerror < greedyError)
                    newstate ~= State(newerror, newvec);
                }
            }
            // writeln(states);

            // newstateの上位だけ取り出してstatesを置き換える
            states.length = min(newstate.length, _selectM);
            topNCopy!"a.error < b.error"(newstate, states);
            newstate.length = 0;
        }

        if(states.length == 0) {
            return GreedyDetector!C().detect(recvY, chMatR, decoded);
        } else {
            // states.sort!"a.error < b.error"();
            auto top = states.reduce!"a.error < b.error ? a : b"();
            foreach(i; 0 .. chMatR.length!1)
                decoded[i] = top.vec[i];

            return top.error;
        }
    }


  private:
    alias F = typeof(C.init.re);

    const(C)[] _points;
    size_t _selectM;

    static struct State
    {
        F error;
        Vector!(C, Contiguous) vec;
    }

}


struct GreedyDetector(C)
{
    this(const(C)[] points)
    {
        _points = points;
    }


    F detect(Vector!(C, Contiguous) recvY, Matrix!(C, Contiguous) chMatR, C[] decoded)
    {
        auto vec = decoded.sliced.vectored;
        F error = 0;
        foreach_reverse(row; 0 .. chMatR.length!0)
        {
            size_t minIndex = 0;
            F minError = F.infinity;

            foreach(i, p; _points) {
                vec[row] = p;
                F err = error + sqAbs(recvY[row] - dot(chMatR.rowVec(row), vec));
                if(err < minError) {
                    minError = err;
                    minIndex = i;
                }
            }

            vec[row] = _points[minIndex];
            error = minError;
        }

        return error;
    }


  private:
    alias F = typeof(C.init.re);

    const(C)[] _points;
}


final class QRDMLDDetector(Algo, C) : IDetector!(C, C)
if(isComplex!C && (is(typeof(C.init.re) == float) || is(typeof(C.init.re) == double)) )
{
    alias InputElementType = C;
    alias OutputElementType = C;

    import std.stdio;

    this(Mat)(Algo algo, in Mat chMat)
    in(isMatrixLike!Mat && chMat.length!0 == chMat.length!1)
    {
        _algo = algo;

        _chMatQ = matrix!C(chMat.length!0, chMat.length!0);
        _chMatR = matrix!C(chMat.length!0, chMat.length!1);
        _recvY = vector!C(chMat.length!0);
        qrDecomp(chMat, _chMatQ, _chMatR);
        _chMatQ[] = _chMatQ.H;
    }


    size_t inputLength() const { return _chMatR.length!0; }
    size_t outputLength() const { return _chMatR.length!1; }


    C[] detect(in C[] received, return ref C[] detected)
    in(received.length % this.inputLength == 0)
    {
        immutable M = this.inputLength;
        immutable N = this.outputLength;

        if(detected.length != received.length / M * N)
            detected.length = received.length / M * N;

        foreach(n; 0 .. received.length / M) {
            _recvY[] = received[n * M .. (n+1) * M].sliced.vectored;
            _recvY[] = _chMatQ * _recvY;
            _algo.detect(_recvY, _chMatR, detected[n * N .. (n+1) * N]);
        }

        return detected;
    }


  private:
    alias F = typeof(C.init.re);

    Algo _algo;
    Vector!(C, Contiguous) _recvY;
    Matrix!(C, Contiguous) _chMatQ;
    Matrix!(C, Contiguous) _chMatR;
}


auto makeQRMMLDDetector(C, Mod, Mat)(Mod mod, in Mat chMat, size_t selectM)
in(mod.symOutputLength == 1)
{
    auto points = mod.allConstellationPoints();
    auto algo = MAlgorithmDetector!C(points, selectM);

    return new QRDMLDDetector!(typeof(algo), C)(algo, chMat);
}


auto makeSphereDetector(C, Mod, Mat, F)(Mod mod, in Mat chMat, F radius = typeof(C.init.re).infinity, Flag!"dynamicUpdate" dynamicUpdate = Yes.dynamicUpdate)
in(mod.symOutputLength == 1)
{
    auto points = mod.allConstellationPoints();
    auto algo = SphereDetector!C(points, radius, dynamicUpdate);

    return new QRDMLDDetector!(typeof(algo), C)(algo, chMat);
}


auto makeQRDGreedyDetector(C, Mod, Mat)(Mod mod, in Mat chMat)
in(mod.symOutputLength == 1)
{
    auto points = mod.allConstellationPoints();
    auto algo = GreedyDetector!C(points);

    return new QRDMLDDetector!(typeof(algo), C)(algo, chMat);
}


unittest
{
    import std.stdio;
    import mir.ndslice : sliced;
    import dffdd.math.exprtemplate;
    import std.math;
    import dffdd.mod.qpsk;

    alias C = MirComplex!float;
    auto chMat = [C(1, 0), C(0, 0), C(0, 0), C(1, 0)].sliced(2, 2).matrixed;
    auto recv = [C(1/SQRT2, 1/SQRT2), C(-1/SQRT2, 1/SQRT2)];

    auto detector1 = makeQRMMLDDetector!C(QPSK!C(), chMat, 20);
    C[] dst1;
    dst1 = detector1.detect(recv, dst1);
    assert(dst1[0].re.isClose(recv[0].re));
    assert(dst1[0].im.isClose(recv[0].im));
    assert(dst1[1].re.isClose(recv[1].re));
    assert(dst1[1].im.isClose(recv[1].im));

    auto detector2 = makeSphereDetector!C(QPSK!C(), chMat);
    C[] dst2;
    dst2 = detector2.detect(recv, dst2);
    assert(dst2[0].re.isClose(recv[0].re));
    assert(dst2[0].im.isClose(recv[0].im));
    assert(dst2[1].re.isClose(recv[1].re));
    assert(dst2[1].im.isClose(recv[1].im));

    auto detector3 = makeQRDGreedyDetector!C(QPSK!C(), chMat);
    C[] dst3;
    dst3 = detector3.detect(recv, dst3);
    assert(dst3[0].re.isClose(recv[0].re));
    assert(dst3[0].im.isClose(recv[0].im));
    assert(dst3[1].re.isClose(recv[1].re));
    assert(dst3[1].im.isClose(recv[1].im));
}
