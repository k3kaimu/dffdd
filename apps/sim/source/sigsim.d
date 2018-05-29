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
    {
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
            txpas = _tempbuf[5][0 .. len];

        // 自端末の送信機のIQインバランス
        if(_model.useSTXIQ) _txIQMixer(xs, txiqs);
        else                txiqs[] = xs[];

        // 自端末の送信機のPAの歪み
        _txPAVGA(txiqs, txvgas);
        if(_model.useSTXPA) _txPARapp(txvgas, txpas);
        else                txpas[] = txvgas[];

        C[] rxants = _tempbuf[6][0 .. len],
            rxvgas = _tempbuf[7][0 .. len],
            rxlnas = _tempbuf[8][0 .. len],
            rxiqs = _tempbuf[9][0 .. len],
            rxqzvgas = _tempbuf[10][0 .. len],
            rxqzs = _tempbuf[11][0 .. len];

        _channel(txpas, rxants);
        _rxLNAVGA(rxants, rxvgas);

        if(_nowTrainingMode){
            foreach(i; 0 .. len)
                rxvgas[i] = rxvgas[i] * _selfInterferenceCoef + ns[i] * _noiseCoef;
        }else{
            foreach(i; 0 .. len)
                rxvgas[i] = rxvgas[i] * _selfInterferenceCoef + ds[i] * _desiredCoef + ns[i] * _noiseCoef;
        }

        // 自端末の受信機のLNAの歪み
        if(_model.useSRXLN) _rxLNARapp(rxvgas, rxlnas);
        else                rxlnas[] = rxvgas[];

        // 自端末の受信機のIQインバランス
        if(_model.useSRXIQ) _rxIQMixer(rxlnas, rxiqs);
        else                rxiqs[] = rxlnas[];

        _rxQZVGA(rxiqs, rxqzvgas);

        // 自端末の受信機のAD変換器の量子化誤差
        if(_model.useSRXQZ) _rxQZ(rxqzvgas, rxqzs);
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

        dst._txIQMixer = this._txIQMixer.dup;
        dst._txPAVGA = this._txPAVGA.dup;
        dst._txPARapp = this._txPARapp.dup;
        
        dst._channel = this._channel.dup;

        dst._rxLNAVGA = this._rxLNAVGA.dup;
        dst._rxLNARapp = this._rxLNARapp.dup;
        dst._rxIQMixer = this._rxIQMixer.dup;
        dst._rxQZVGA = this._rxQZVGA.dup;
        dst._rxQZ = this._rxQZ.dup;

        return dst;
    }


    JSONValue info()
    {
        JSONValue dst = JSONValue(string[string].init);
        dst["txIQMixer"] = _txIQMixer.dumpInfoToJSON();
        dst["txPAVGA"] = _txPAVGA.dumpInfoToJSON();
        dst["txPARapp"] = _txPARapp.dumpInfoToJSON();
        dst["channel"] = _channel.dumpInfoToJSON();
        dst["rxLNAVGA"] = _rxLNAVGA.dumpInfoToJSON();
        dst["rxLNARapp"] = _rxLNARapp.dumpInfoToJSON();
        dst["rxIQMixer"] = _rxIQMixer.dumpInfoToJSON();
        dst["rxQZVGA"] = _rxQZVGA.dumpInfoToJSON();
        dst["rxQZ"] = _rxQZ.dumpInfoToJSON();
        return dst;
    }


    bool isConvergedAllVGA() @property
    {
        static
        bool _isconverged(X)(X* v)
        {
            static if(is(typeof((X x){ return x.isConverged; })))
                return v.isConverged;
            else
                return true;
        }

        bool dst = true;

        // .isConvergedをメンバとして持つメンバ変数を列挙してチェック
        foreach(m; __traits(allMembers, typeof(this)))
            static if(is(typeof(mixin(m))) && is(typeof(&mixin(m)) == typeof(mixin(m))*))
                dst = dst && _isconverged(&mixin(m));

        return dst;
    }


    C[] linearSIChannel() @property
    {
        C[] channel = _channel.coefficients.dup;

        Gain g = Gain.fromPowerGain(1);
        g *= _txIQMixer.gain;
        g *= _txPAVGA.gain;
        g *= _txPARapp.linearGain;
        g *= _rxLNAVGA.gain;
        g *= _rxLNARapp.linearGain;
        g *= _rxIQMixer.gain;
        g *= _rxQZVGA.gain;

        foreach(ref e; channel)
            e *= g.gain;

        return channel;
    }


  private:
    Model _model;

  public:
    Signal desiredBaseband;
    Signal txBaseband;
    Signal noise;

  private:
    C[][13] _tempbuf;

    bool _nowTrainingMode;
    bool* _useSWPOFDM;

    C _desiredCoef = C(1);
    C _selfInterferenceCoef = C(1);
    C _noiseCoef = C(1);

    IQImbalanceConverter!C _txIQMixer;
    PowerControlAmplifierConverter!C _txPAVGA;
    RappModelConverter!C _txPARapp;

    FIRFilterConverter!C _channel;

    PowerControlAmplifierConverter!C _rxLNAVGA;
    RappModelConverter!C _rxLNARapp;
    IQImbalanceConverter!C _rxIQMixer;
    PowerControlAmplifierConverter!C _rxQZVGA;
    SimpleQuantizerConverter!C _rxQZ;
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
                        .map!"a*1.0L".asSignal;

    dst.txBaseband = siRandomBits(model)
                        .connectToModulator(modOFDM(model), dst._useSWPOFDM, model)
                        .map!"a*1.0L".asSignal;

    dst.noise = thermalNoise(model).connectTo!VGA(model.lna.NF).toWrappedRange;

    dst._txIQMixer = IQImbalanceConverter!C(0.dB, model.txIQMixer.imbCoef);
    dst._txPAVGA = PowerControlAmplifierConverter!C((model.pa.TX_POWER.dBm - model.pa.GAIN.dB).dBm, 1e-2);
    dst._txPARapp = RappModelConverter!C(model.pa.GAIN, 1, model.pa.IIP3.volt / 2);

    dst._channel = FIRFilterConverter!C(model.channel.impulseResponse[0 .. model.channel.taps]);
    dst._rxLNAVGA = PowerControlAmplifierConverter!C(model.thermalNoise.power(model) * model.lna.NF * model.INR, 1e-2);
    dst._rxLNARapp = RappModelConverter!C(model.lna.GAIN, model.lna.smoothFactor, (model.lna.IIP3.dBm - 36).dB.gain);
    dst._rxIQMixer = IQImbalanceConverter!C(0.dB, model.rxIQMixer.imbCoef);
    dst._rxQZVGA = PowerControlAmplifierConverter!C((30 - model.ofdm.PAPR.dB + 4.76).dBm, 1e-2);
    dst._rxQZ = SimpleQuantizerConverter!C(model.quantizer.numOfBits);

    return dst;
}
