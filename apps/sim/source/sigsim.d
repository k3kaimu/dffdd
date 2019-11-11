module sigsim;

import std.algorithm;
import std.complex;
import std.file;
import std.format;
import std.math;
import std.meta;
import std.path;
import std.range;
import std.traits;
import std.typecons;
import std.json;

import dffdd.mod.qam;
import dffdd.blockdiagram.adder;
import dffdd.blockdiagram.amplifier;
import dffdd.blockdiagram.utils;
import dffdd.blockdiagram.iqmixer;
import dffdd.blockdiagram.filter;
import dffdd.blockdiagram.quantizer;
import dffdd.utils.unit;
import dffdd.utils.json;

import models;
import simmain;
import snippet;


alias Signal = ForwardRange!(Complex!real);


Signal asSignal(R)(R r)
if(isForwardRange!R)
{
    return r.toWrappedRange;
}


final class SimulatedSignals
{
    alias C = Complex!real;


    this()
    {
        foreach(ref buf; _tempbuf)
            buf = new C[1024];
    }


    void popFrontN(size_t n)
    {
        assert(_tempbuf[0].length != 0);

        if(n <= _tempbuf[0].length) {
            this.fillBuffer!(["txBaseband"])(_tempbuf[$-1][0 .. n]);
        }else{
            while(n != 0){
                size_t consumed = min(n, _tempbuf[0].length);
                this.popFrontN(consumed);
                n -= consumed;
            }
        }
    }


    enum bool empty = false;


    void trainAGC()
    out {
        import std.meta : AliasSeq;
        foreach(obj; AliasSeq!(_txPAVGA, _rxLNAVGA, _rxQZVGA))
            if(!obj.isNull)
                assert(obj.get.isConverged);
    }
    do {
        _nowTrainingMode = true;

        immutable oldUseSWP = *_useSWPOFDM;
        *_useSWPOFDM = false;

        while(!this.isConvergedAllVGA){
            this.popFrontN(_model.ofdm.numOfSamplesOf1Symbol);
        }

        _nowTrainingMode = false;
        *_useSWPOFDM = oldUseSWP;

        this.popFrontN(_model.ofdm.numOfSamplesOf1Symbol * 4);
    }


    void fillBuffer(alias ms, X, size_t N)(X[][N] buffers...)
    if(ms.length == N && N != 0)
    in{
        auto s = buffers[0].length;
        foreach(e; buffers) assert(e.length == s);
    }
    body{
        static foreach(m; aliasSeqOf!ms)
            static assert(canFind(["txBaseband", "desiredBaseband", "paDirect", "received"], m), m ~ " is illegal.");

        immutable size_t len = buffers[0].length;
        if(len == 0) return;

        foreach(ref buf; _tempbuf)
            buf.length = max(buf.length, len);

        C[] ds = _tempbuf[0][0 .. len],
            xs = _tempbuf[1][0 .. len],
            ns = _tempbuf[2][0 .. len];

        foreach(i; 0 .. len){
            ds[i] = desiredBaseband.front;
            xs[i] = txBaseband.front;
            ns[i] = noise.front;

            desiredBaseband.popFront();
            txBaseband.popFront();
            noise.popFront();
        }

        C[] txiqs = _tempbuf[3][0 .. len],
            txvgas = _tempbuf[4][0 .. len],
            txpas = _tempbuf[5][0 .. len],
            detxiqs = _tempbuf[6][0 .. len],
            detxvgas = _tempbuf[7][0 .. len],
            detxpas = _tempbuf[8][0 .. len];

        // 自端末の送信機のIQインバランス
        if(_model.useSTXIQ) _txIQMixer.get()(xs, txiqs);
        else                txiqs[] = xs[];

        // 相手端末の送信機のIQインバランス
        if(_model.useDTXIQ) _detxIQMixer.get()(ds, detxiqs);
        else                detxiqs[] = ds[];

        // 自端末の送信機のPAの歪み
        _txPAVGA.get()(txiqs, txvgas);
        if(_model.useSTXPA) _txPANonlin.get()(txvgas, txpas);
        else                txpas[] = txvgas[];

        // 相手端末の送信機のPAの歪み
        _detxPAVGA.get()(detxiqs, detxvgas);
        if(_model.useDTXPA) _detxPANonlin.get()(detxvgas, detxpas);
        else                detxpas[] = detxvgas[];

        C[] rxants = _tempbuf[9][0 .. len],
            rxvgas = _tempbuf[10][0 .. len],
            rxlnas = _tempbuf[11][0 .. len],
            rxiqs = _tempbuf[12][0 .. len],
            rxqzvgas = _tempbuf[13][0 .. len],
            rxqzs = _tempbuf[14][0 .. len],
            rxds = _tempbuf[15][0 .. len],
            rxdsants = _tempbuf[16][0 .. len];

        _channelSI.get()(txpas, rxants);
        _rxLNAVGA.get()(rxants, rxvgas);
        _channelDesired.get()(detxpas, rxdsants);
        _rxDESVGA.get()(rxdsants, rxds);

        if(_nowTrainingMode){
            foreach(i; 0 .. len)
                rxvgas[i] = rxvgas[i] * _selfInterferenceCoef + ns[i] * _noiseCoef;
        }else{
            foreach(i; 0 .. len)
                rxvgas[i] = rxvgas[i] * _selfInterferenceCoef + rxds[i] * _desiredCoef + ns[i] * _noiseCoef;
        }

        // 自端末の受信機のLNAの歪み
        if(_model.useSRXLN) _rxLNANonlin.get()(rxvgas, rxlnas);
        else                rxlnas[] = rxvgas[];

        // 自端末の受信機のIQインバランス
        if(_model.useSRXIQ) _rxIQMixer.get()(rxlnas, rxiqs);
        else                rxiqs[] = rxlnas[];

        _rxQZVGA.get()(rxiqs, rxqzvgas);

        // 自端末の受信機のAD変換器の量子化誤差
        if(_model.useSRXQZ) _rxQZ.get()(rxqzvgas, rxqzs);
        else                rxqzs[] = rxqzvgas[];

        static foreach(i, m; aliasSeqOf!ms) {
            static if(m == "txBaseband")
                xs.copy(buffers[i]);
            else static if(m == "desiredBaseband")
                ds.copy(buffers[i]);
            else static if(m == "paDirect")
                txpas.copy(buffers[i]);
            else static if(m == "received")
                rxqzs.copy(buffers[i]);
            else static assert(0, m ~ " is illegal.");
        }
    }


    void ignoreDesired(bool b) @property
    {
        _desiredCoef = b ? C(0) : C(1);
    }


    void ignoreSI(bool b) @property
    {
        _selfInterferenceCoef = b ? C(0) : C(1);
    }


    void ignoreNoise(bool b) @property
    {
        _noiseCoef = b ? C(0) : C(1);
    }


    void useSWPOFDM(bool b) @property
    {
        *_useSWPOFDM = b;
    }


    SimulatedSignals save() @property
    {
        SimulatedSignals dst = new SimulatedSignals();
        dst._model = this._model;
        
        dst.desiredBaseband = this.desiredBaseband.save;
        dst.txBaseband = this.txBaseband.save;
        dst.noise = this.noise.save;

        foreach(i, e; _tempbuf)
            dst._tempbuf[i].length = e.length;

        dst._nowTrainingMode = this._nowTrainingMode;
        dst._useSWPOFDM = new bool(*this._useSWPOFDM);

        dst._desiredCoef = this._desiredCoef;
        dst._selfInterferenceCoef = this._selfInterferenceCoef;
        dst._noiseCoef = this._noiseCoef;

        if(!this._txIQMixer.isNull)     dst._txIQMixer = this._txIQMixer.get.dup;
        if(!this._txPAVGA.isNull)       dst._txPAVGA = this._txPAVGA.get.dup;
        if(!this._txPANonlin.isNull)    dst._txPANonlin = this._txPANonlin.get.dup;

        if(!this._detxIQMixer.isNull)     dst._detxIQMixer = this._detxIQMixer.get.dup;
        if(!this._detxPAVGA.isNull)       dst._detxPAVGA = this._detxPAVGA.get.dup;
        if(!this._detxPANonlin.isNull)    dst._detxPANonlin = this._detxPANonlin.get.dup;
        
        if(!this._channelSI.isNull)         dst._channelSI = this._channelSI.get.dup;
        if(!this._channelDesired.isNull)    dst._channelDesired = this._channelDesired.get.dup;

        if(!this._rxDESVGA.isNull)      dst._rxDESVGA = this._rxDESVGA.get.dup;
        if(!this._rxLNAVGA.isNull)      dst._rxLNAVGA = this._rxLNAVGA.get.dup;
        if(!this._rxLNANonlin.isNull)   dst._rxLNANonlin = this._rxLNANonlin.get.dup;
        if(!this._rxIQMixer.isNull)     dst._rxIQMixer = this._rxIQMixer.get.dup;
        if(!this._rxQZVGA.isNull)       dst._rxQZVGA = this._rxQZVGA.get.dup;
        if(!this._rxQZ.isNull)          dst._rxQZ = this._rxQZ.get.dup;

        return dst;
    }


    JSONValue info()
    {
        JSONValue dst = JSONValue(string[string].init);
        if(!_txIQMixer.isNull)  dst["txIQMixer"] = _txIQMixer.get.dumpInfoToJSON();
        if(!_txPAVGA.isNull)    dst["txPAVGA"] = _txPAVGA.get.dumpInfoToJSON();
        if(!_txPANonlin.isNull)   dst["txPANonlin"] = _txPANonlin.get.dumpInfoToJSON();
        if(!_detxIQMixer.isNull)  dst["detxIQMixer"] = _detxIQMixer.get.dumpInfoToJSON();
        if(!_detxPAVGA.isNull)    dst["detxPAVGA"] = _detxPAVGA.get.dumpInfoToJSON();
        if(!_detxPANonlin.isNull)   dst["detxPANonlin"] = _detxPANonlin.get.dumpInfoToJSON();
        if(!_channelSI.isNull)    dst["channelSI"] = _channelSI.get.dumpInfoToJSON();
        if(!_channelDesired.isNull)    dst["channelDesired"] = _channelDesired.get.dumpInfoToJSON();
        if(!_rxLNAVGA.isNull)   dst["rxLNAVGA"] = _rxLNAVGA.get.dumpInfoToJSON();
        if(!_rxLNANonlin.isNull)  dst["rxLNARapp"] = _rxLNANonlin.get.dumpInfoToJSON();
        if(!_rxIQMixer.isNull)  dst["rxIQMixer"] = _rxIQMixer.get.dumpInfoToJSON();
        if(!_rxQZVGA.isNull)    dst["rxQZVGA"] = _rxQZVGA.get.dumpInfoToJSON();
        if(!_rxQZ.isNull)       dst["rxQZ"] = _rxQZ.get.dumpInfoToJSON();
        return dst;
    }


    bool isConvergedAllVGA() @property
    {
        static
        bool _isconverged(X)(X* v)
        {
            static if(is(X : Nullable!Y, Y)){
                if(v.isNull)
                    return true;
                else
                    return _isconverged(&(v.get()));
            }else{
                static if(is(typeof((X x){ return x.isConverged; })))
                    return v.isConverged;
                else
                    return true;
            }
        }

        bool dst = true;

        // .isConvergedをメンバとして持つメンバ変数を列挙してチェック
        foreach(m; __traits(allMembers, typeof(this)))
            static if(is(typeof(mixin(m))) && is(typeof(&mixin(m)) == typeof(mixin(m))*))
                dst = dst && _isconverged(&mixin(m));

        return dst;
    }


    // Gain desiredSignalGain() @property
    // {
    //     Gain g = Gain.fromPowerGain(1);
    //     if(!_rxDESVGA.isNull)       g *= _rxDESVGA.gain;
    //     if(!_rxLNANonlin.isNull)    g *= _rxLNANonlin.linearGain;
    //     if(!_rxIQMixer.isNull)      g *= _rxIQMixer.gain;
    //     if(!_rxQZVGA.isNull)        g *= _rxQZVGA.gain;
    //     return g;
    // }


    C[] linearSIChannel() @property
    {
        C[] channel = _channelSI.get.coefficients.dup;

        Gain g = Gain.fromPowerGain(1);
        if(!_txIQMixer.isNull)  g *= _txIQMixer.get.gain;
        if(!_txPAVGA.isNull)    g *= _txPAVGA.get.gain;
        if(!_txPANonlin.isNull)   g *= _txPANonlin.get.linearGain;
        if(!_rxLNAVGA.isNull)   g *= _rxLNAVGA.get.gain;
        if(!_rxLNANonlin.isNull)  g *= _rxLNANonlin.get.linearGain;
        if(!_rxIQMixer.isNull)  g *= _rxIQMixer.get.gain;
        if(!_rxQZVGA.isNull)    g *= _rxQZVGA.get.gain;

        foreach(ref e; channel)
            e *= g.asV;

        return channel;
    }


  private:
    Model _model;

  public:
    Signal desiredBaseband;
    Signal txBaseband;
    Signal noise;

  private:
    C[][17] _tempbuf;

    bool _nowTrainingMode;
    bool* _useSWPOFDM;

    C _desiredCoef = C(1);
    C _selfInterferenceCoef = C(1);
    C _noiseCoef = C(1);

    Nullable!(IQImbalanceConverter!C) _txIQMixer;
    Nullable!(PowerControlAmplifierConverter!C) _txPAVGA;
    Nullable!(RappModelConverter!C) _txPANonlin;

    Nullable!(IQImbalanceConverter!C) _detxIQMixer;
    Nullable!(PowerControlAmplifierConverter!C) _detxPAVGA;
    Nullable!(RappModelConverter!C) _detxPANonlin;

    Nullable!(FIRFilterConverter!C) _channelSI;
    Nullable!(FIRFilterConverter!C) _channelDesired;

    Nullable!(PowerControlAmplifierConverter!C) _rxDESVGA;
    Nullable!(PowerControlAmplifierConverter!C) _rxLNAVGA;
    Nullable!(RappModelConverter!C) _rxLNANonlin;
    Nullable!(IQImbalanceConverter!C) _rxIQMixer;
    Nullable!(PowerControlAmplifierConverter!C) _rxQZVGA;
    Nullable!(SimpleQuantizerConverter!C) _rxQZ;
}


SimulatedSignals makeSimulatedSignals(Model model, string resultDir = null)
{
    alias C = Complex!real;

    immutable doOutput = resultDir !is null;
    immutable bool* alwaysFalsePointer = new bool(false);
    immutable bool* alwaysTruePointer = new bool(true);

    SimulatedSignals dst = new SimulatedSignals();
    dst._model = model;
    dst._useSWPOFDM = new bool(false);

    dst.desiredBaseband = desiredRandomBits(model)
                        .connectToModulator(modOFDM(model), alwaysFalsePointer, model)
                        .map!"a*1.0L"
                        .connectTo!(PowerControlAmplifierConverter!C)(30.dBm, 1e-2)
                        .asSignal;

    dst.txBaseband = siRandomBits(model)
                        .connectToModulator(modOFDM(model), dst._useSWPOFDM, model)
                        .map!"a*1.0L"
                        .connectTo!(PowerControlAmplifierConverter!C)(30.dBm, 1e-2)
                        .asSignal;

    dst.noise = thermalNoise(model).connectTo!VGA(model.lna.NF).toWrappedRange;

    if(model.useSTXIQ)
        dst._txIQMixer = IQImbalanceConverter!C(0.dB, model.txIQMixer.imbCoef);

    if(model.useDTXIQ)
        dst._detxIQMixer = IQImbalanceConverter!C(0.dB, model.txIQMixer.imbCoef);

    if(model.useSTXPA){
        dst._txPAVGA = PowerControlAmplifierConverter!C(model.pa.TX_POWER / model.pa.GAIN, 1e-2);
        dst._txPANonlin = RappModelConverter!C(model.pa.GAIN, model.pa.smoothFactor, model.pa.IIP3.volt / 2);
    }else{
        dst._txPAVGA = PowerControlAmplifierConverter!C(model.pa.TX_POWER, 1e-2);
    }

    if(model.useDTXPA) {
        dst._detxPAVGA = PowerControlAmplifierConverter!C(model.pa.TX_POWER / model.pa.GAIN, 1e-2);
        dst._detxPANonlin = RappModelConverter!C(model.pa.GAIN, model.pa.smoothFactor, model.pa.IIP3.volt / 2);
    }else{
        dst._detxPAVGA = PowerControlAmplifierConverter!C(model.pa.TX_POWER, 1e-2);
    }

    dst._channelSI = FIRFilterConverter!C(model.channelSI.impulseResponse[0 .. model.channelSI.taps]);
    dst._channelDesired = FIRFilterConverter!C(model.channelDesired.impulseResponse[0 .. model.channelDesired.taps]);
    
    Voltage receivedSIPower = model.thermalNoise.power(model) * model.lna.NF * model.INR;
    if(model.useSRXLN){
        dst._rxLNAVGA = PowerControlAmplifierConverter!C(receivedSIPower, 1e-2);
        dst._rxLNANonlin = RappModelConverter!C(model.lna.GAIN, model.lna.smoothFactor, (model.lna.IIP3 / 36.dBm).asV);
        //  dst._rxLNANonlin = SoftLimitConverter!C(model.lna.GAIN, (model.lna.IIP3 / 36.dBm).asV);
        // dst._rxLNANonlin = SalehModelConverter!C(model.lna.GAIN, (model.lna.IIP3 / 36.dBm).asV);
    }else{
        dst._rxLNAVGA = PowerControlAmplifierConverter!C(receivedSIPower, 1e-2);
        // dst._rxLNANonlin = RappModelConverter!C(model.lna.GAIN, model.lna.smoothFactor, real.infinity);
    }

    if(model.useSRXIQ)
        dst._rxIQMixer = IQImbalanceConverter!C(0.dB, model.rxIQMixer.imbCoef);

    dst._rxQZVGA = PowerControlAmplifierConverter!C((30 - model.ofdm.PAPR.asdB + 4.76).dBm, 1e-2);

    if(model.useSRXQZ)
        dst._rxQZ = SimpleQuantizerConverter!C(model.quantizer.numOfBits);

    // immutable snScaleOFDM = Gain.fromPowerGain(1.0L * model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling / model.ofdm.numOfSubcarrier);
    // immutable noisePowerPerSubcarrier = model.thermalNoise.power(model) * model.lna.NF / snScaleOFDM;

    dst._rxDESVGA = PowerControlAmplifierConverter!C(model.thermalNoise.power(model) * model.lna.NF * model.SNR, 1e-2);

    return dst;
}
