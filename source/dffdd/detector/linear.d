module dffdd.detector.linear;

import mir.ndslice : slice, sliced, Contiguous;
import dffdd.math.matrix;
import dffdd.math.linalg;
import dffdd.math.vector;
import dffdd.math.matrixspecial;

import dffdd.mod.primitives;


final class LinearDetector(C, Mod)
{
    alias InputElementType = C;
    alias OutputElementType = Bit;


    this(M)(Mod mod, M mat)
    if(isMatrixLike!M)
    {
        _mod = mod;
        _mat = matrix!(M.ElementType)(mat.length!0, mat.length!1);
        _mat[] = mat;
        _filted = new M.ElementType[](mat.length!0);
    }


    size_t inputLength() const { return _mat.length!1; }
    size_t outputLength() const { return _mat.length!0 * _mod.symInputLength * _mod.symOutputLength; }


    Bit[] detect(in C[] received, return ref Bit[] detected)
    in(received.length % this.inputLength == 0)
    {
        immutable M = this.outputLength;
        immutable N = this.inputLength;

        if(detected.length != received.length / N * M)
            detected.length = received.length / N * M;

        foreach(i; 0 .. received.length / this.inputLength) {
            auto recvY = received[i * N .. (i+1) * N].sliced.vectored;
            _filted.sliced.vectored()[] = _mat * recvY;

            auto s = detected[i * M .. (i+1) * M];
            _mod.demodulate(_filted, s);
        }

        return detected;
    }


  private:
    Mod _mod;
    C[] _filted;
    Matrix!(C, Contiguous) _mat;
}


auto makeMMSEDetector(C, Mod, Mat)(Mod mod, Mat chMat, double sigma2)
{
    immutable M = chMat.length!0;
    immutable N = chMat.length!1;

    auto invMat = matrix!(Mat.ElementType)(M, M);
    invMat[] = sigma2 * identity!(Mat.ElementType)(M) + chMat * chMat.H;
    invMat[] = inv(invMat);
    
    auto wmat = matrix!(Mat.ElementType)(N, M);
    wmat[] = chMat.H * invMat;

    return new LinearDetector!(C, Mod)(mod, wmat);
}


auto makeZFDetector(C, Mod, Mat)(Mod mod, Mat chMat)
{
    immutable M = chMat.length!0;
    immutable N = chMat.length!1;

    // auto invMat = matrix!C(Mat.ElementType)(chMat)
    auto invMat = matrix!(Mat.ElementType)(N, N);
    invMat[] = chMat.H * chMat;
    invMat[] = inv(invMat);

    auto wmat = matrix!(Mat.ElementType)(N, M);
    wmat[] = invMat * chMat.H;

    return new LinearDetector!(C, Mod)(mod, wmat);
}

