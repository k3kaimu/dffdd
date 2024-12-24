/++ dub.json: {
    "name": "sigre",
    "authors": [
        "Kazuki Komatsu"
    ],
    "copyright": "Copyright © 2024, Kazuki Komatsu",
    "description": "Signal Regenerator",
    "license": "NYSL",
    "dependencies": {
        "dffdd": { "path": "../.." }
    },
    "lflags-linux": ["-L/opt/OpenBLAS/lib"],
    "libs": ["openblas", "lapacke"],
}
+/


/// BFList X[3,0] L[1,0] M[1,2,3]
/// X[p,q] := x^p conj(x[n])^q
/// L[p,q] := laguerreOBF(x[n], p, q)
/// M[a,*b,c,d,*e,...] := x[n-a] * conj(x[n-b]) * x[n-c] * x[n-d] * conj(x[n-e]) * ...

/// $ sigre <nTrain> <nRegenBlockSize> <Regenerator Type> [Regenerator Parameters...]
/// $ sigre <nTrain> <nRegenBlockSize> <PHammLS/CHammLS> <nTaps> @ <Basis Function Type> [Basis Functions Parameters...]
/// $ sigre 1000 2000 PHammLS 32 @ BFList X[1,0] X[0,1] X[3,0] X[2,1] X[1,2] X[3,0] # 32タップで(x, x^{*}, x^3, x|x|^2, x^{*}|x|^2, x^{*}^3)を基底にもつキャンセラ
/// $ sigre 1000 2000 PHammLS 32 @ BFList L[3,0] L[2,1] L[1,2] # 32タップで(L_{3,0}, L_{2,1}, L_{1,2})を基底にもつキャンセラ
/// $ sigre 1000 2000 PHammLS 32 @ BFList M[] M[0,1], M[0,*1] M[0,*0,1] # 32タップでメモリをもつ4つの基底（1, x[n] * x[n-1], x[n] * x^{*}[n-1], x[n] * x^{*}[n] * x[n-1]）からなるキャンセラ
/// $ sigre 1000 2000 PHammLS 32 @ <Basis Function Type> [Basis Functions Paramters] % Orthonormalizer # 後続する基底関数を直交化する
/// $ sigre 1000 2000 PHammLS 32 @ OFDMUpSampler 64 16 8 % <Basis Function Type> [Basis Functions Paramters] % OFDMDownSampler 64 16 8 # 入力信号がFFTサイズ64，CP長16のOFDM信号だと仮定して，それを8倍にオーバーサンプリングしてから基底関数にかけて，1/8倍にデシメーションして元に戻してからキャンセラに入れる
import std;
import dffdd;
import std.complex : conj, sqAbs;
import mir.ndslice : slice, sliced;


immutable string paramSeparator = "@";
immutable string blockSeparator = "%";


immutable string helpText =
`
OVERVIEW: SigRe - Signal Regenerator

USAGE: sigre <nTrain> <nRegenBlockSize> <Regenerator Type> [Regenerator Parameters...] @ <Basis Function Builders...>

Examples:
$ sigre 1000 80 PHammLS 32 @ OFDMUpSampler 64 16 8 % BFList X[1,0] L[0,1] X[3,0] M[0,1,*1] % OFDMDownSampler 64 16 8 % Orthonormalizer
`;


void main(string[] args)
{
    if(args.length == 1) {
        writeln(helpText);
        return;
    }

    mainImpl(args[1 .. $], stdin, stdout);
}


void mainImpl(string[] args, File input, File output)
{
    alias C = Complex!float;

    string[] params = args;
    immutable nTrain = params[0].to!ulong.failureMsg!"a > 0"("The number of training samples must be larger than zero.");
    immutable nRegenBlockSize = params[1].to!ulong.failureMsg!"a > 0"("The number of regenerating block samples must be larger than zero.");

    C[] xsTrain = new C[](nTrain);
    C[] ysTrain = new C[](nTrain);
    enforce(input.rawRead(xsTrain).length == nTrain, "Cannot read input samples for training.");
    enforce(input.rawRead(ysTrain).length == nTrain, "Cannot read output samples for training.");

    // C[] xsRegen = new C[](nRegenBlockSize);
    // enforce(input.rawRead(xsRegen).length == nTrain, "Cannot read input samples for regeneration.");

    immutable regenName = params[2];
    immutable cntRegenParams = params[3 .. $].countUntil("@");
    string[] regenParams  = params[3 .. 3 + cntRegenParams];
    string[] basisParams = params[4 + cntRegenParams .. $];

    IRegenerator!C regenImpl;
    switch(regenName) {
        case "PHammLS":
            regenImpl = new PHammLSRegenerator!C(xsTrain, ysTrain, regenParams, makeBBB!C(basisParams));
            break;
        case "CHammLS":
            break;
        default:
            enforce(0, "Undefined regenerator type: '%s'".format(regenName));
    }

    C[] xsRegen = new C[](nRegenBlockSize);
    while(! input.eof) {
        size_t readed = input.rawRead(xsRegen).length;
        if(readed != nRegenBlockSize) {
            enforce(readed == 0, "insufficient input signal");
            break;
        }

        output.rawWrite(regenImpl.regenerate(xsRegen));
    }
}


version(SIGRE_TEST_MODE) static this()
{
    alias C = Complex!float;

    size_t nSamples = 1000;
    size_t nBlockSize = 200;

    C[] xs = [];
    foreach(i; 0 .. nSamples) {
        xs ~= C(uniform01(), uniform01());
    }

    // C[] hs = [C(-1, 0.5), C(0.1, 2), C(3, 0), C(0, -1)];
    C[] hs = [C(-1, 0.5), C(0.1, 2)];

    C[] ys = new C[nSamples];
    foreach(i; 0 .. nSamples) {
        ys[i] = C(0, 0);
        foreach(k; 0 .. hs.length) {
            ys[i] += hs[k] * (i >= k ? xs[i - k] : C(0));
            ys[i] += hs[k] * (i >= k ? xs[i - k]^^2 * conj(xs[i - k]) : C(0));
        }
    }

    scope(exit) { std.file.remove("input"); std.file.remove("output"); }

    {
        File infile = File("input", "w");
        infile.rawWrite(xs);
        infile.rawWrite(ys);
        infile.rawWrite(xs);
    }
    
    mainImpl([nSamples.to!string, nBlockSize.to!string, "PHammLS", "4", "@", "Orthonormalizer", "%", "BFList", "X[1,0]", "X[2,1]", "%", "Orthonormalizer"], File("input", "r"), File("output", "w"));

    C[] regen = File("output", "r").rawRead(new C[nSamples]);
    foreach(k; 0 .. nSamples / nBlockSize) {
        auto buff = regen[k * nBlockSize + 5 .. (k+1) * nBlockSize];
        auto ysff = ys[k * nBlockSize + 5 .. (k+1) * nBlockSize];
        assert(zip(buff, ysff).map!"a[0] - a[1]".map!sqAbs.sum / (nBlockSize - 5) < 1e-6);
    }

    assert(0);
}


interface IRegenerator(C)
{
    C[] regenerate(in C[] xs);
}


final class PHammLSRegenerator(C) : IRegenerator!C
{
    this(in C[] trs, in C[] ysTrain, string[] params, IBasisBuildingBlock!C bbb)
    {
        _bbb = bbb;
        C[][] xsTrain = _bbb.generate([trs]);

        _nBF = xsTrain.length;

        enforce(xsTrain.length > 0, "[PHammLS] The dimension of input signal for training must be larger than 1.");
        enforce(params.length > 0, "[PHammLS] PHammLS needs nTaps as an argument.");
        _nTaps = params[0].to!uint.failureMsg("[PHammLS] invalid number of taps '%s'".format(params[0]));

        immutable nTrainSamples = xsTrain[0].length;
        enforce(nTrainSamples > _nTaps, "[PHammLS] The number of training samples must be larger than nTaps.");
        enforce(xsTrain[0].length == ysTrain.length, "[PHammLS] The number of transmission and reception samples for training must be same.");

        _estCoeffs = estimateCoeffs(xsTrain, ysTrain, params);
    }


    C[] regenerate(in C[] input)
    {
        auto xs = _bbb.generate([input]);
        immutable nRegenSamples = xs[0].length;
        auto ysRegen = slice!C(nRegenSamples).vectored;
        auto xsMatrix = slice!C(nRegenSamples, _nBF * _nTaps).matrixed;

        foreach(i; 0 .. nRegenSamples) {
            foreach(k; 0 .. _nBF) {
                foreach(j; 0 .. _nTaps) {
                    xsMatrix[i, k * _nTaps + j] = (i >= j) ? xs[k][i - j] : C(0, 0);
                }
            }
        }

        ysRegen[] = xsMatrix * _estCoeffs.sliced.vectored;

        C[] dst = new C[ysRegen.length];
        foreach(i; 0 .. ysRegen.length)
            dst[i] = ysRegen[i];

        return dst;
    }


  private:
    size_t _nBF, _nTaps;
    C[] _estCoeffs;
    IBasisBuildingBlock!C _bbb;

    C[] estimateCoeffs(in C[][] xsTrain, in C[] ysTrain, /*in C[][] xsRegen,*/ string[] params)
    {
        immutable nTrainSamples = xsTrain[0].length;
        auto mx = slice!C(nTrainSamples - _nTaps, _nBF * _nTaps);
        foreach(i; _nTaps .. nTrainSamples) {
            foreach(k; 0 .. _nBF) {
                foreach(j; 0 .. _nTaps) {
                    mx[i - _nTaps, k * _nTaps + j] = xsTrain[k][i - j];
                }
            }
        }

        C[] ysT = new C[nTrainSamples - _nTaps];
        ysT[] = ysTrain[_nTaps .. $];
        auto est = leastSquareEstimateRowMajor(mx, ysT);
        return est;
    }
}


// unittest
version(SIGRE_TEST_MODE) static this()
{
    alias C = Complex!float;

    C[] xs = [C(1, 0), C(2, 1), C(1, 0), C(3, 3), C(-1, -5)];
    C[] ys = [C(1, 0)*2, C(2, 1)*2, C(1, 0)*2, C(3, 3)*2, C(-1, -5)*2];
    C[] xs2 = [C(1, 1), C(2, 0), C(-1, 5), C(3, 0), C(-1, -1)];

    auto regen = new PHammLSRegenerator!C(xs, ys, 5, ["1"], new NullBBB!C());

    C[] ys2 = regen.regenerate(xs2);
    foreach(i; 0 .. xs2.length) {
        assert(approxEqualC(xs2[i] * 2, ys2[i]));
    }
}

version(SIGRE_TEST_MODE) static this()
{
    alias C = Complex!float;

    C[] xs;
    foreach(i; 0 .. 1000) {
        xs ~= C(uniform01(), uniform01());
    }

    C[][] hs = [
        [C(-1, 0.5), C(0.1, 2), C(3, 0), C(0, -1)],
        [C(0, 0.2), C(0.4, 3), C(2, 2), C(1, -1)],
    ];

    C[] ys = new C[1000];
    foreach(i; 0 .. 1000) {
        ys[i] = C(0, 0);
        foreach(k; 0 .. 4) {
            ys[i] += hs[0][k] * (i >= k ? xs[i - k] : C(0));
            ys[i] += hs[1][k] * (i >= k ? xs[i - k]^^2 : C(0));
        }
    }

    auto regenImpl = new PHammLSRegenerator!C(xs, ys, 5, ["4"], makeBBB!C(["BFList", "X[1,0]", "X[2,0]"]));

    C[] regen = regenImpl.regenerate(xs);
    foreach(i; 0 .. regen.length) {
        assert(approxEqualC(regen[i], ys[i]));
    }
}


template failureMsg(conds...)
{
    T failureMsg(T)(lazy T value, string msg)
    {
        T dst;
        try {
            dst = value();
        } catch(Exception ex) {
            enforce(false, msg ~ '\n' ~ ex.to!string);
        }

        foreach(c; conds)
            enforce(unaryFun!c(dst), msg);

        return dst;
    }
}


interface IBasisBuildingBlock(C)
{
    C[][] generate(in C[][] xs);
}


final class BFList(C) : IBasisBuildingBlock!C
{
    this(string[] params)
    {
        _parsed = parseBasisIndexList(params);
        _originalParams = params;
    }


    C[][] generate(in C[][] xs)
    {
        enforce(xs.length == 1, "[BFList] Invalid input dimension which must be 1.");
        return generateImpl(xs[0]);
    }


    C[][] generateImpl(C)(in C[] xs)
    {
        immutable nSamples = xs.length;

        C[][] applied = new C[][](_parsed.length, nSamples);
        foreach(i, p; _parsed) {
            switch(p.label) {
                case 'X':
                    enforce(p.value.length == 2, "[BFList.X] '%s' is invalid: The size of the index-pair must be two.".format(_originalParams[i]));
                    enforce(p.options == [null, null], "[BFList.X] '%s' is invalid: The options are ignored.".format(_originalParams[i]));
                    foreach(j; 0 .. nSamples)
                        applied[i][j] = xs[j]^^(p.value[0]) * conj(xs[j])^^(p.value[1]);
                    
                    break;

                case 'L':
                    enforce(p.value.length == 2, "[BFList.L] '%s' is invalid: The size of the index-pair must be two.".format(_originalParams[i]));
                    enforce(p.options == [null, null], "[BFList.L] '%s' is invalid: The options are ignored.".format(_originalParams[i]));
                    foreach(j; 0 .. nSamples)
                        applied[i][j] = laguerreOBF(xs[j], p.value[0], p.value[1]);

                    break;

                case 'M':
                    foreach(j; 0 .. nSamples) {
                        applied[i][j] = C(1, 0);
                        foreach(k; 0 .. p.value.length) {
                            enforce(p.options[k] == "" || p.options[k] == "*", "[BFList.M] '%s' is invalid. Each option must be '' or '*'.".format(_originalParams[i]));

                            immutable C v = (j >= p.value[k]) ? xs[j - p.value[k]] : C(0, 0);
                            applied[i][j] *= p.options[k] == "*" ? conj(v) : v;
                        }
                    }
                    break;

                default:
                    enforce(0, "[BFList] Invalid basis function '%s'.".format(_originalParams[i]));
            }
        }

        return applied;
    }

  private:
    IndexTuple[] _parsed;
    string[] _originalParams;

    static struct IndexTuple
    {
        char label;
        int[] value;
        string[] options;
    }

    static IndexTuple[] parseBasisIndexList(string[] params)
    {
        IndexTuple[] ret;
        foreach(i, strpair; params) {
            if(strpair == blockSeparator || strpair == paramSeparator) break;

            enforce(std.ascii.isAlpha(strpair[0]) && strpair[1] == '[' && strpair[$-1] == ']', "[Basis] The %s-th index pair '%s' is invalid format.".format(i, strpair));

            auto ss = strpair[2 .. $-1].split(",");
            IndexTuple parsed = IndexTuple(strpair[0], new int[](ss.length), new string[](ss.length));

            foreach(k; 0 .. ss.length) {
                if(ss[k][0] == '*') {
                    parsed.options[k] = "*";
                    ss[k] = ss[k][1 .. $];
                } else {
                    parsed.options[k] = null;
                }

                auto n = ss[k].to!int.failureMsg("[Basis] The %s-th index pair '%s' is invalid format. The valid format is (p, q,...) for integers p, g, ....".format(i, strpair));
                parsed.value[k] = n;
            }

            ret ~= parsed;
        }

        params = params[ret.length .. $];
        return ret;
    }
}


// unittest
version(SIGRE_TEST_MODE) static this()
{
    alias C = Complex!float;

    string[] args = ["X[1,0]", "X[0,1]", "X[3,0]", "L[4,2]", "M[0,1,*1,2,2,*2]", "M[]"];
    C[] xs = [C(1, 0), C(2, 1), C(3, 2)];

    C[][] dst = new BFList!C(args).generate([xs]);

    assert(dst[0][0] == xs[0]);
    assert(dst[0][1] == xs[1]);
    assert(dst[0][2] == xs[2]);

    assert(dst[1][0] == xs[0].conj);
    assert(dst[1][1] == xs[1].conj);
    assert(dst[1][2] == xs[2].conj);

    assert(dst[2][0] == xs[0]^^3);
    assert(dst[2][1] == xs[1]^^3);
    assert(dst[2][2] == xs[2]^^3);

    assert(dst[3][0] == laguerreOBF(xs[0], 4, 2));
    assert(dst[3][1] == laguerreOBF(xs[1], 4, 2));
    assert(dst[3][2] == laguerreOBF(xs[2], 4, 2));

    assert(dst[4][0] == 0);
    assert(dst[4][1] == 0);
    assert(dst[4][2] == xs[2] * xs[1] * conj(xs[1]) * xs[0] * xs[0] * conj(xs[0]));

    assert(dst[5][0] == C(1, 0));
    assert(dst[5][1] == C(1, 0));
    assert(dst[5][2] == C(1, 0));
}


final class Orthonormalizer(C) : IBasisBuildingBlock!C
{
    this() {}

    C[][] generate(in C[][] xs)
    {
        bool isInit = false;
        if(_coefs is null) {
            _coefs = new C[][](xs.length, xs.length);
            _norms = new C[](xs.length);
            isInit = true;
        }

        return gramschmidt(xs, isInit);
    }

  private:
    C[][] _coefs;
    C[] _norms;

    C[][] gramschmidt(in C[][] xs, bool isInit)
    {
        if(isInit) {
            _norms[] = C(1);
            foreach(i, ref e; _coefs) {
                e[] = C(0);
                e[i] = C(1);
            }
        }

        immutable nBF = xs.length;
        immutable N = xs[0].length;

        C[][] dst = new C[][](nBF, N);

        static C ip(in C[] a, in C[] b) {
            C sum = C(0, 0);
            foreach(i; 0 .. a.length) sum += conj(a[i]) * b[i];
            return sum / a.length;
        }

        foreach(i; 0 .. nBF) {
            dst[i][] = xs[i][];
            foreach(k; 0 .. i) {
                if(isInit)  _coefs[i][k] = ip(dst[k], dst[i]);
                immutable c = _coefs[i][k];
                foreach(n; 0 .. N)
                    dst[i][n] -= c * dst[k][n];
            }

            if(isInit) _norms[i] = sqrt(ip(dst[i], dst[i]).re);
            immutable norm = _norms[i];
            foreach(n; 0 .. N)
                dst[i][n] /= norm;
        }

        return dst;
    }
}


// unittest
version(SIGRE_TEST_MODE) static this()
{
    alias C = Complex!float;

    C[][] xs = [
        [C(1, 0), C(2, 1), C(3, 2)],
        [C(2, 0), C(1, 1), C(8, 7)]
    ];

    C[][] dst = new Orthonormalizer!C().generate(xs);
    
    C sum = 0;
    foreach(i; 0 .. xs[0].length)
        sum += dst[0][i] * conj(dst[1][i]);

    assert(sum.re < 1e-6);
    assert(sum.im < 1e-6);

    foreach(i; 0 .. dst.length) {
        double power = 0;
        foreach(e; dst[i])
            power += std.complex.sqAbs(e);

        assert((power / dst[i].length).approxEqual(1));
    }
}


final class OFDMUpSampler(C : Cpx!F, alias Cpx, F) : IBasisBuildingBlock!C
{
    this(uint nFFT, uint nCP, uint nOS)
    {
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;

        _fftL = makeFFTWObject!Cpx(nFFT);
        _fftH = makeFFTWObject!Cpx(nFFT * nOS);
    }


    C[][] generate(in C[][] xs)
    {
        C[][] dst = new C[][](xs.length);
        foreach(i; 0 .. xs.length)
            dst[i] = genImpl(xs[i]);

        return dst;
    }

  private:
    size_t _nFFT, _nCP, _nOS;
    FFTWObject!Cpx _fftH, _fftL;

    C[] genImpl(in C[] xs)
    {
        enforce(xs.length % (_nFFT + _nCP) == 0, "[OFDMUpSampler] The length of the input signal (%s) must be a constant multiple of _nFFT + _nCP = %s.".format(xs.length, _nFFT + _nCP));
    
        immutable nSamples = xs.length;
        immutable nSym = xs.length / (_nFFT + _nCP);

        C[] dst = new C[]((_nFFT + _nCP) * _nOS * nSym);
        foreach(i; 0 .. xs.length / (_nFFT + _nCP)) {
            _fftL.inputs!F[] = xs[(_nFFT + _nCP) * i + _nCP .. (_nFFT + _nCP) * (i+1)];
            _fftL.fft!F();
            _fftH.inputs!F[] = C(0, 0);
            _fftH.inputs!F[0 .. _nFFT / 2] = _fftL.outputs!F[0 .. _nFFT / 2];
            _fftH.inputs!F[$ - _nFFT/2 .. $] = _fftL.outputs!F[_nFFT / 2 .. $];
            _fftH.ifft!F();

            auto dstsym = dst[i * (_nFFT + _nCP) * _nOS * nSym .. (i+1) * (_nFFT + _nCP) * _nOS * nSym];
            dstsym[0 .. _nCP * _nOS] = _fftH.outputs!F[$ - _nCP * _nOS .. $];
            dstsym[_nCP * _nOS .. $] = _fftH.outputs!F[];
        }

        // 元の電力に戻す
        immutable scale = sqrt(_nOS * 1.0);
        foreach(ref e; dst) e *= scale;

        return dst;
    }
}


final class OFDMDownSampler(C : Cpx!F, alias Cpx, F) : IBasisBuildingBlock!C
{
    this(uint nFFT, uint nCP, uint nOS)
    {
        _nFFT = nFFT;
        _nCP = nCP;
        _nOS = nOS;

        _fftL = makeFFTWObject!Cpx(nFFT);
        _fftH = makeFFTWObject!Cpx(nFFT * nOS);
    }


    C[][] generate(in C[][] xs)
    {
        C[][] dst = new C[][](xs.length);
        foreach(i; 0 .. xs.length)
            dst[i] = genImpl(xs[i]);

        return dst;
    }

  private:
    size_t _nFFT, _nCP, _nOS;
    FFTWObject!Cpx _fftH, _fftL;


    C[] genImpl(C)(in C[] xs)
    {
        enforce(xs.length % ((_nFFT + _nCP)*_nOS) == 0, "[OFDMDownSampler] The length of the input signal (%s) must be a constant multiple of (_nFFT + _nCP)*_nOS = %s.".format(xs.length, (_nFFT + _nCP)*_nOS));
        
        immutable nSamples = xs.length;
        immutable nSym = xs.length / ((_nFFT + _nCP) * _nOS);

        C[] dst = new C[]((_nFFT + _nCP) * nSym);
        foreach(i; 0 .. xs.length / ((_nFFT + _nCP)*_nOS) ) {
            _fftH.inputs!F[] = xs[(_nFFT + _nCP) * _nOS * i + _nCP * _nOS .. (_nFFT + _nCP) * _nOS * (i+1)];
            _fftH.fft!F();
            _fftL.inputs!F[0 .. _nFFT / 2] = _fftH.outputs!F[0 .. _nFFT / 2];
            _fftL.inputs!F[_nFFT/2 .. $] = _fftH.outputs!F[$ - _nFFT / 2 .. $];
            _fftL.ifft!F();

            auto dstsym = dst[i * (_nFFT + _nCP) * nSym .. (i+1) * (_nFFT + _nCP) * nSym];
            dstsym[0 .. _nCP] = _fftL.outputs!F[$ - _nCP .. $];
            dstsym[_nCP .. $] = _fftL.outputs!F[];
        }

        // 元の電力に戻す
        immutable scale = sqrt(1.0 / _nOS);
        foreach(ref e; dst) e *= scale;

        return dst;
    }
}


// unittest
version(SIGRE_TEST_MODE) static this()
{
    alias C = Complex!float;

    C[] xs = [C(1, 0), C(2, 1), C(1, 0)];

    C[][] ys = new OFDMDownSampler!C(2, 1, 8).generate(new OFDMUpSampler!C(2, 1, 8).generate([xs]));
    foreach(i; 0 .. xs.length) assert(approxEqualC(xs[i], ys[0][i]));
}


final class ChainedBBB(C) : IBasisBuildingBlock!C
{
    this(IBasisBuildingBlock!C pre, IBasisBuildingBlock!C post)
    {
        _pre = pre;
        _post = post;
    }


    C[][] generate(in C[][] xs)
    {
        return _post.generate(_pre.generate(xs));
    }


  private:
    IBasisBuildingBlock!C _pre;
    IBasisBuildingBlock!C _post;
}


final class NullBBB(C) : IBasisBuildingBlock!C
{
    this() {}


    C[][] generate(in C[][] xs)
    {
        C[][] dst = new C[][](xs.length, xs[0].length);
        foreach(i; 0 .. xs.length)
            dst[i][] = xs[i][];
        
        return dst;
    }
}


IBasisBuildingBlock!C makeBBB(C)(string[] params)
{
    if(params.length == 0) return new NullBBB!C();

    Tuple!(string, "name", string[], "params") getCurrentBBBNameAndParams()
    {
        if(params[0] == blockSeparator) params.popFront();

        string name = params[0];
        params.popFront();

        string[] ret;
        while(params.length != 0 && params[0] != blockSeparator) {
            ret ~= params[0];
            params.popFront();
        }

        return typeof(return)(name, ret);
    }

    IBasisBuildingBlock!C ret = new NullBBB!C();
    while(params.length != 0) {
        auto curr = getCurrentBBBNameAndParams();

        switch(curr.name) {
            case "BFList":
                ret = new ChainedBBB!C(ret, new BFList!C(curr.params));
                break;                

            case "Orthonormalizer":
                enforce(curr.params.length == 0, "[Orthonormalizer] Given parameters %s are ignored.".format(curr.params));
                ret = new ChainedBBB!C(ret, new Orthonormalizer!C());
                break;

            case "OFDMUpSampler":
                enforce(curr.params.length == 3, "[OFDMUpSampler] Given parameters %s are ignored.".format(curr.params[3 .. $]));
                immutable nFFT = curr.params[0].to!uint.failureMsg!"a > 0"("[OFDMUpSampler] This function needs the parameter of FFT points as a positive integer.");
                immutable nCP = curr.params[1].to!uint.failureMsg!"a >= 0"("[OFDMUpSampler] This function needs the parameter of the CP points as a positive integer.");
                immutable nOS = curr.params[2].to!uint.failureMsg!"a > 0"("[OFDMUpSampler] This function needs the parameter of the oversampling factor as a positive integer.");
                ret = new ChainedBBB!C(ret, new OFDMUpSampler!C(nFFT, nCP, nOS));
                break;

            case "OFDMDownSampler":
                enforce(curr.params.length == 3, "[OFDMDownSampler] Given parameters %s are ignored.".format(curr.params[3 .. $]));
                immutable nFFT = curr.params[0].to!uint.failureMsg!"a > 0"("[OFDMDownSampler] This function needs the parameter of FFT points as a positive integer.");
                immutable nCP = curr.params[1].to!uint.failureMsg!"a >= 0"("[OFDMDownSampler] This function needs the parameter of the CP points as a positive integer.");
                immutable nOS = curr.params[2].to!uint.failureMsg!"a > 0"("[OFDMDownSampler] This function needs the parameter of the oversampling factor as a positive integer.");
                ret = new ChainedBBB!C(ret, new OFDMDownSampler!C(nFFT, nCP, nOS));
                break;

            default:
                enforce(false, "Invalid BasisFunction or Modifier name '%s'".format(params[0]));
                return null;
        }
    }

    return ret;
}
