module dffdd.detector.tridiag_dp;

import dffdd.math.complex;
import dffdd.math.linalg;
import dffdd.math.matrix;
import dffdd.math.vector;

import mir.ndslice : Contiguous, sliced;

import dffdd.mod.primitives;
import dffdd.detector.primitives;


/** 
 * エルミート行列からなる通信路行列を3重対角化し，動的計画法により最尤推定します
 * ただし，送信ベクトルは三重対角に使用するユニタリ行列でプリコーディング及びポストコーディングされている必要があります．
 */
final class TridiagDPMLDetector(C) : IDetector!(C, C)
if(isComplex!C && (is(typeof(C.init.re) == float) || is(typeof(C.init.re) == double)) )
{
    alias InputElementType = C;
    alias OutputElementType = C;

    import std.stdio;

    this(Tridiagonalized!C tdres, const(C)[] points)
    {
        _tdres = tdres;
        _points = points;
        _recvY = vector!C(tdres.q.length!0);
        _states = new State[][](_points.length, _points.length);
        _newstates = new State[][](_points.length, _points.length);

        foreach(i; 0 .. _points.length) foreach(j; 0 .. _points.length) {
            _states[i][j] = State(0, vector!C(_tdres.q.length!0));
            _newstates[i][j] = State(0, vector!C(_tdres.q.length!0));
        }
    }


    size_t inputLength() const { return _tdres.q.length!0; }
    size_t outputLength() const { return _tdres.q.length!0; }


    C[] detect(in C[] received, return ref C[] detected)
    in(received.length % this.inputLength == 0)
    {
        immutable M = this.inputLength;
        immutable N = this.outputLength;

        if(detected.length != received.length / M * N)
            detected.length = received.length / M * N;

        foreach(n; 0 .. received.length / M) {
            _recvY[] = received[n * M .. (n+1) * M].sliced.vectored;
            // _recvY[] = _tdres.q.H * _recvY;
            detectImpl(detected[n * N .. (n+1) * N]);
        }

        return detected;
    }


  private:
    alias F = typeof(C.init.re);

    const(C)[] _points;
    Tridiagonalized!C _tdres;
    Vector!(C, Contiguous) _recvY;
    State[][] _states, _newstates;


    void detectImpl(C[] detected)
    {
        import std.algorithm : swap;

        immutable N = detected.length;

        // 最初の行の計算
        foreach(i; 0 .. _points.length) foreach(j; 0 .. _points.length) {
            _newstates[i][j].err = F.infinity;
            _newstates[i][j].vec[] = C(0);
            _states[i][j].err = 0;
            _states[i][j].vec[] = C(0);
            _states[i][j].vec[0] = _points[i];
            _states[i][j].vec[1] = _points[j];
            _states[i][j].err = sqAbs(_recvY[0] - _tdres.diag[0] * _points[i] - _tdres.offdiag[0] * _points[j]);
        }


        // 最初と最後の行以外の行についての計算
        foreach(r; 1 .. N - 1) {
            foreach(i; 0 .. _points.length) foreach(j; 0 .. _points.length) {
                _newstates[i][j].err = F.infinity;
                _newstates[i][j].vec[] = C(0);
            }

            foreach(i; 0 .. _points.length) foreach(j; 0 .. _points.length) foreach(k; 0 .. _points.length) {
                F newerr = _states[i][j].err;
                newerr += sqAbs(_recvY[r] - _tdres.offdiag[r-1] * _points[i]
                                         - _tdres.diag[r] * _points[j]
                                        - _tdres.offdiag[r] * _points[k]);
                
                if(newerr < _newstates[j][k].err) {
                    _newstates[j][k].err = newerr;
                    _newstates[j][k].vec[] = _states[i][j].vec;
                    _newstates[j][k].vec[r+1] = _points[k];
                }
            }

            swap(_newstates, _states);
        }

        // 最終行の計算
        State minState;
        minState.err = F.infinity;

        foreach(i; 0 .. _points.length) foreach(j; 0 .. _points.length) {
            _states[i][j].err += sqAbs(_recvY[N-1] - _tdres.offdiag[N-2] * _points[i] - _tdres.diag[N-1] * _points[j]);

            if(_states[i][j].err < minState.err) {
                minState = _states[i][j];
            }
        }

        foreach(i; 0 .. N)
            detected[i] = minState.vec[i];
    }


    static struct State
    {
        F err;
        Vector!(C, Contiguous) vec;
    }
}


auto makeTridiagDPMLDetector(C, Mod)(Mod mod, Tridiagonalized!C tdres)
in(mod.symOutputLength == 1)
{
    auto points = mod.allConstellationPoints();

    return new TridiagDPMLDetector!(C)(tdres, points);
}



unittest
{
    import std.stdio;
    import mir.ndslice : sliced;
    import dffdd.math.exprtemplate;
    import std.math;
    import dffdd.mod.qpsk;
    import dffdd.math.matrixspecial;

    alias C = MirComplex!float;
    // auto chMat = [C(1, 0), C(0, 0), C(0, 0), C(1, 0)].sliced(2, 2).matrixed;

    auto chMat = matrix!C(4, 4);
    chMat[] = identity!C(4);
    auto recv = [C(1/SQRT2, 1/SQRT2), C(1/SQRT2, -1/SQRT2), C(-1/SQRT2, 1/SQRT2), C(-1/SQRT2, -1/SQRT2)];

    auto tdres = tridiagonalize(chMat);

    auto detector1 = makeTridiagDPMLDetector!C(QPSK!C(), tdres);
    C[] dst1;
    dst1 = detector1.detect(recv, dst1);

    foreach(i; 0 .. 4) {
        assert(dst1[i].re.isClose(recv[i].re));
        assert(dst1[i].im.isClose(recv[i].im));
    }
}
