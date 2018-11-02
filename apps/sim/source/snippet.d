module snippet;


import std.stdio;
import std.math;
import std.complex;
import std.range;
import std.algorithm;
import std.traits;

import dffdd.dsp.statistics;
import carbon.math;

void writePSD(R)(auto ref R r, File file, real samplingFreq, size_t resolution)
{
    auto psd = r.calculatePowerSpectralDensity(samplingFreq, resolution);
    foreach(i, e; psd){
        file.writefln("%f,%f,%f,", (i*1.0/resolution*samplingFreq-(samplingFreq/2))/1e6, e, e == 0 ? -400 : 30+10*log10(e));
    }
}


real measurePower(R)(auto ref R r, size_t n)
{
    real s = 0;
    foreach(i; 0 .. n){
        if(r.empty){
            n = i;
            break;
        }

        auto e = r.front;
        s += e.re^^2 + e.im^^2;
        r.popFront();
    }

    return s / n;
}


real measureBER(R1, R2)(ref R1 r1, ref R2 r2, ulong totalBits)
{
    r1.popFrontN(10*1000);
    r2.popFrontN(10*1000);

    ulong cnt;
    foreach(i; 0 .. totalBits){
        auto a = r1.front;
        auto b = r2.front;

        if(a != b) ++cnt;

        r1.popFront();
        r2.popFront();
    }

    return (cast(real)cnt)/totalBits;
}


auto connectToSwitch(R)(R r, const(bool)* sw, ElementType!R zero = complexZero!(ElementType!R))
{
    static struct Result 
    {
        auto front() @property
        {
            if(*_sw)
                return _r.front;
            else
                return _zero;
        }


        void popFront()
        {
            _r.popFront();
        }


        bool empty()
        {
            return _r.empty;
        }


      static if(isForwardRange!R)
      {
        typeof(this) save() @property
        {
            typeof(return) dst = this;

            dst._r = this._r.save;

            return dst;
        }
      }


      private:
        R _r;
        const(bool)* _sw;
        ElementType!R _zero;
    }


    return Result(r, sw, zero);
}


final class OFDMSymbolExchanger(Mod)
{
    alias InputElementType = Mod.InputElementType;
    alias OutputElementType = Mod.OutputElementType;


    this(Mod mod, const(bool)* sw, uint nFFT, uint nCp, uint nTone, uint nUpSampling = 1)
    {
        _mod = mod;
        _sw = sw;
        _nFFT = nFFT;
        _nCp = nCp;
        _nTone = nTone;
        _nUpSampling = nUpSampling;
    }


    size_t symInputLength() const @property { return _mod.symInputLength*2; }
    size_t symOutputLength() const @property { return _mod.symOutputLength*2; }


    ref OutputElementType[] modulate(in InputElementType[] input, return ref OutputElementType[] output)
    in{
        assert(input.length % this.symInputLength == 0);
    }
    body{
        _mod.modulate(input, output);

        if(*_sw){
            foreach(i; 0 .. input.length / this.symInputLength) {
                auto symAB = output[i * this.symOutputLength .. (i+1)*this.symOutputLength];

                auto startOfA = symAB[_nCp * _nUpSampling .. $];
                auto startOfB = symAB[$/2 + _nCp * _nUpSampling .. $];

                foreach(j; 0 .. _nFFT * _nUpSampling / 2)
                    swap(startOfA[j], startOfB[j]);
            }
        }

        return output;
    }


    ref InputElementType[] demodulate(in OutputElementType[] input, return ref InputElementType[] output)
    {
        if(output.length != input.length / this.symOutputLength * this.symInputLength)
            output.length = input.length / this.symOutputLength * this.symInputLength;

        auto input2 = input.dup;

        if(*_sw){
            foreach(i; 0 .. input.length / this.symOutputLength) {
                auto symAB = input2[i * this.symOutputLength .. (i+1)*this.symOutputLength];

                auto startOfA = symAB[_nCp * _nUpSampling .. $];
                auto startOfB = symAB[$/2 + _nCp * _nUpSampling .. $];

                foreach(j; 0 .. _nFFT * _nUpSampling / 2)
                    swap(startOfA[j], startOfB[j]);
            }
        }

        _mod.demodulate(input2, output);

        return output;
    }


  private:
    Mod _mod;
    const(bool)* _sw;
    uint _nFFT, _nCp, _nTone, _nUpSampling;
}


OFDMSymbolExchanger!Mod makeOFDMSymbolExchanger(Mod)(Mod mod, const(bool)* sw, uint nFFT, uint nCp, uint nTone, uint nUpSampling = 1)
{
    return new OFDMSymbolExchanger!Mod(mod, sw, nFFT, nCp, nTone, nUpSampling);
}


auto loggingTo(R, O)(R r, O orange)
// if(isForwardRange!R && isOutputRange!(O, R))
{
    static struct Result
    {
        this(R r, O o, bool isOwner)
        {
            _r = r;
            _o = o;
            _isOwner = isOwner;
            _isOutputted = false;
        }


        ElementType!R front()
        {
            auto f = _r.front;
            if(!_isOutputted && _isOwner){
                .put(_o, f);
                _isOutputted = true;
            }

            return f;
        }


        void popFront()
        {
            _r.popFront();
            _isOutputted = false;
        }


        bool empty() { return _r.empty; }


        Result save()
        {
            typeof(return) dst = this;

            dst._isOwner = false;
            dst._r = this._r.save;

            return dst;
        }


      private:
        bool _isOwner;
        bool _isOutputted;
        R _r;
        O _o;
    }


    return Result(r, orange, true);
}


template checkConcept(alias concept)
{
    auto checkConcept(R)(R r)
    if(concept!R)
    {
        return r;
    }
}



static
final class ModulatedRange(R, Mod, bool isDemod = false)
// if(isForwardRange!R)
{
  static if(isDemod)
  {
      alias A = Mod.OutputElementType;
      alias B = Mod.InputElementType;

      private size_t inputLen() { return _mod.symOutputLength; }
      private size_t outputLen() { return _mod.symInputLength; }
  }
  else
  {
      alias A = Mod.InputElementType;
      alias B = Mod.OutputElementType;

      private size_t inputLen() { return _mod.symInputLength; }
      private size_t outputLen() { return _mod.symOutputLength; }
  }


    private this() {}


    this(R r, Mod modObj)
    {
        _r = r;
        _mod = modObj;
        _inputBuf = new A[this.inputLen];
        _outputBuf = new B[this.outputLen];
        _idx = this.outputLen - 1;
        this.popFront();
    }


    B front() @property { return _outputBuf[_idx]; }
    bool empty() @property { return _idx >= this.outputLen; }


    void popFront()
    {
        if(_idx < this.outputLen - 1){
            ++_idx;
            return;
        }

        _idx = 0;
        foreach(i, ref e; _inputBuf){
            if(_r.empty) goto Lempty;
            e = _r.front;
            _r.popFront();
        }

      static if(isDemod)
        _mod.demodulate(_inputBuf, _outputBuf);
      else
        _mod.modulate(_inputBuf, _outputBuf);

        return;

      Lempty:
        _idx = this.outputLen;
    }


  static if(isForwardRange!R)
  {
    typeof(this) save() @property
    {
        typeof(return) dst = new ModulatedRange();

        dst._r = this._r.save;
        dst._mod = this._mod;
        dst._inputBuf = this._inputBuf.dup;
        dst._outputBuf = this._outputBuf.dup;
        dst._idx = this._idx;

        return dst;
    }
  }

    private:
    R _r;
    Mod _mod;
    A[] _inputBuf;
    B[] _outputBuf;
    size_t _idx;
}


final class NopCanceller(C)
{
    this() {}

    enum inputBlockLength = 1;

    void apply(Flag!"learning" doLearning)(in C[] input, in C[] received, C[] errors)
    in{
        assert(input.length == received.length);
        assert(input.length == errors.length);
    }
    do{
        errors[] = received[];
    }


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator) {}
}


auto makeSignalLoggingCanceller(C, Canceller)(Canceller canceller)
{
    static
    final class SignalLoggingImpl
    {
        import std.base64;
        import std.json;

        enum inputBlockLength = 1;

        this(Canceller canc)
        {
            _canc = canc;
        }

        void apply(Flag!"learning" doLearning)(in C[] input, in C[] received, C[] errors)
        in{
            assert(input.length == received.length);
            assert(input.length == errors.length);
        }
        body{
            // errors[] = received[];
            _canc.apply!doLearning(input, received, errors);

            if(doLearning){
                _ssTrX ~= input;
                _ssTrY ~= received;
                _ssTrZ ~= errors;
            }else{
                _ssCnX ~= input;
                _ssCnY ~= received;
                _ssCnZ ~= errors;
            }
        }


        JSONValue info() @property
        {
            JSONValue[string] dst;
            dst["ssTrX_Base64"] = Base64.encode(cast(ubyte[])_ssTrX);
            dst["ssTrY_Base64"] = Base64.encode(cast(ubyte[])_ssTrY);
            dst["ssTrZ_Base64"] = Base64.encode(cast(ubyte[])_ssTrZ);
            dst["ssCnX_Base64"] = Base64.encode(cast(ubyte[])_ssCnX);
            dst["ssCnY_Base64"] = Base64.encode(cast(ubyte[])_ssCnY);
            dst["ssCnZ_Base64"] = Base64.encode(cast(ubyte[])_ssCnZ);
            return JSONValue(dst);
        }


        void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator) {}


      private:
        Canceller _canc;
        // 学習用信号と，除去用信号
        C[] _ssTrX, _ssTrY, _ssTrZ; // X = 送信, Y = 受信, Z = 除去後
        C[] _ssCnX, _ssCnY, _ssCnZ;
    }

    return new SignalLoggingImpl(canceller);
}



/*
final class SignalLoggingCanceller(C)
{
    import std.base64;
    import std.json;


    this() {}

    enum inputBlockLength = 1;


    void apply(Flag!"learning" doLearning)(in C[] input, in C[] received, C[] errors)
    in{
        assert(input.length == received.length);
        assert(input.length == errors.length);
    }
    body{
        errors[] = received[];

        if(doLearning){
            _ssTrX ~= input;
            _ssTrY ~= received;
        }else{
            _ssCnX ~= input;
            _ssCnY ~= received;
        }
    }


    JSONValue info() @property
    {
        JSONValue[string] dst;
        dst["ssTrX_Base64"] = Base64.encode(cast(ubyte[])_ssTrX);
        dst["ssTrY_Base64"] = Base64.encode(cast(ubyte[])_ssTrY);
        dst["ssCnX_Base64"] = Base64.encode(cast(ubyte[])_ssCnX);
        dst["ssCnY_Base64"] = Base64.encode(cast(ubyte[])_ssCnY);
        return JSONValue(dst);
    }


    void preLearning(M, Signals)(M model, Signals delegate(M) signalGenerator) {}


    // 学習用信号と，除去用信号
    C[] _ssTrX, _ssTrY;
    C[] _ssCnX, _ssCnY;
}
*/