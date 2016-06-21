module dffdd.blockdiagram.mod.ofdm;


import std.string;
import std.conv;

import std.numeric;
import std.complex;
import std.algorithm;

import dranges.range;
import dranges.algorithm;


/*
struct OFDMEqualizer(Mod)
{
  private:
    Mod _mod;
}
*/
struct OFDMEqualizer
{
    static
    auto makeBlock(R, Mod, Rnd)(R r, Mod mod, uint nFFT, uint nCP, uint nSC, uint nOS, Rnd rnd)
    {
        return OFDMEqualizerImpl!(R, Mod, Rnd)(r, mod, nFFT, nCP, nSC, nOS, rnd);
    }


    static struct OFDMEqualizerImpl(R, Mod, Rnd)
    {
        this(R r, Mod mod, uint nFFT, uint nCP, uint nSC, uint nOS, Rnd rnd)
        {
            _r = r;
            _fftObj = new Fft(nFFT * nOS);
            _mod = mod;
            _nFFT = nFFT;
            _nCP = nCP;
            _nSC = nSC;
            _nOS = nOS;
            _rnd = rnd.splitN(_mod.symInputLength);
            _symbol.length = _mod.symOutputLength;
            _buf0.length = _mod.symOutputLength;
            _buf1.length = nFFT * nOS;
            _buf2.length = nFFT * nOS;
            _h.length = nFFT * nOS;
            _empty = false;
        }


        const(cfloat)[] front() const @property
        {
            return _symbol;
        }


        bool empty() const @property { return _empty; }


        void popFront()
        {
            import std.stdio;
            _buf0 = _mod.modulate(_rnd.front, _buf0);
            writeln(_mod.symInputLength);           // 52 * 4
            writeln(_rnd.front.length);             // 
            writeln(_mod.symOutputLength);
            writeln(_buf0.length);
            assert(_buf0.length == _mod.symOutputLength);

            _rnd.popFront();

            _symbol[] = _r.front[];
            _r.popFront();

            foreach(idx; 0 .. _buf0.length / ((_nFFT + _nCP) * _nOS))
            {
                immutable symSize = (_nFFT+_nCP)*_nOS;
                immutable baseIdx = idx*symSize;

                _fftObj.fft(_buf0[baseIdx .. baseIdx + _nFFT * _nOS].map!(a => complex(a.re, a.im)), _buf1);
                _fftObj.fft(_symbol[baseIdx .. baseIdx + _nFFT * _nOS].map!(a => complex(a.re, a.im)), _buf2);

                foreach(i; 0 .. _nFFT * _nOS){
                    //_h[i] = _buf2[i] / (_buf1[i] + 0.001);
                    //_buf2[i] /= _h[i];
                    _buf2[i] = _buf1[i] + 0.001;
                }

                _fftObj.inverseFft(_buf2, _buf1);

                foreach(i, e; _buf1)
                    _symbol[baseIdx + _nCP*_nOS + i] = e.re + e.im*1i;

                _symbol[baseIdx .. baseIdx + _nCP*_nOS] = _symbol[baseIdx + symSize - (_nCP * _nOS) .. baseIdx + symSize];
            }
        }


      private:
        R _r;
        Fft _fftObj;
        Mod _mod;
        uint _nFFT, _nCP, _nSC, _nOS;
        typeof(Rnd.init.splitN(1)) _rnd;
        cfloat[] _symbol;
        cfloat[] _buf0;
        Complex!float[] _buf1, _buf2;
        Complex!float[] _h;
        bool _empty;
    }
}


//import dffdd.utils.msgpackrpc;

/+
struct OFDMModMRPC(string sendName, string recvName)
{
    static
    auto makeBlock(R)(R r, TCPClient client)
    if(isInputRange!R && is(ElementType!R == bool))
    {
        return r.packBinaryDigits!uint(8).connect!(MsgpackPRCBlock!(uint, float[], sendName, recvName, uint.sizeof, cfloat.sizeof))(client).map!"a[0]+a[1]*1i";
    }
}


struct OFDMDemodMRPC(string sendName, string recvName)
{
    static
    auto makeBlock(R)(R r, TCPClient client)
    if(isInputRange!R && is(ElementType!R == bool))
    {
        return r.map!"[a.re, a.im]".connect!(MsgpackPRCBlock!(float[], uint, sendName, recvName, cfloat.sizeof, uint.sizof))(client).toBinaryDigits(8);
    }
}+/


//struct OFDMMod(Mod)