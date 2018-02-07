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

    Tuple!(
        C, "desiredBaseband",
        C, "txBaseband",
        C, "paDirect",
        C, "received"
    ) front() @property
    {
        if(_frontIsComputed) return _cache;

        C d = desiredBaseband.front,
          x = txBaseband.front,
          n = noise.front;

        C txiq, txvga, txpa;
        _txIQMixer(x, txiq);
        _txPAVGA(txiq, txvga);
        _txPARapp(txvga, txpa);

        C rxant, rxvga, rxlna, rxiq, rxqzvga, rxqz;
        _channel(txpa, rxant);
        _rxLNAVGA(rxant, rxvga);

        if(_nowTrainingMode)
            rxvga = rxvga * _selfInterferenceCoef + n * _noiseCoef;
        else
            rxvga = rxvga * _selfInterferenceCoef + d * _desiredCoef + n * _noiseCoef;

        _rxLNARapp(rxvga, rxlna);
        _rxIQMixer(rxlna, rxiq);
        _rxQZVGA(rxiq, rxqzvga);
        _rxQZ(rxqzvga, rxqz);

        typeof(return) dst;
        dst.desiredBaseband = d;
        dst.txBaseband = x;
        dst.paDirect = txpa;
        dst.received = rxqz;
        _cache = dst;
        return dst;
    }


    void popFront()
    {
        if(!_frontIsComputed) this.front();
        desiredBaseband.popFront();
        txBaseband.popFront();
        noise.popFront();
        _frontIsComputed = false;
    }


    enum bool empty = false;


    void trainAGC()
    {
        _nowTrainingMode = true;

        immutable oldUseSWP = *_useSWPOFDM;
        *_useSWPOFDM = false;

        foreach(i; 0 .. _model.numOfModelTrainingSymbols * _model.ofdm.numOfSamplesOf1Symbol)
            this.popFront();

        _nowTrainingMode = false;
        *_useSWPOFDM = oldUseSWP;

        foreach(i; 0 .. _model.numOfModelTrainingSymbols * _model.ofdm.numOfSamplesOf1Symbol)
            this.popFront();
    }


    void fillBuffer(alias ms, X, size_t N)(X[][N] buffers...)
    if(ms.length == N)
    in{
        auto s = buffers[0].length;
        foreach(e; buffers) assert(e.length == s);
    }
    body{
        foreach(i; 0 .. buffers[0].length){
            auto v = this.front;
            foreach(j, m; aliasSeqOf!ms)
                mixin(format(`buffers[j][i] = cast(C)v.%s;`, m));
            
            this.popFront();
        }
    }


    void ignoreDesired(bool b) @property
    {
        _desiredCoef = b ? C(0) : C(1);
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

        dst._cache = this._cache;
        dst._frontIsComputed = this._frontIsComputed;
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


  private:
    Model _model;

  public:
    Signal desiredBaseband;
    Signal txBaseband;
    Signal noise;

  private:
    typeof(front()) _cache;
    bool _frontIsComputed;

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

    SimulatedSignals dst = new SimulatedSignals;
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
    dst._txPAVGA = PowerControlAmplifierConverter!C((model.pa.TX_POWER.dBm - model.pa.GAIN.dB).dBm);
    dst._txPARapp = RappModelConverter!C(model.pa.GAIN, 1, model.pa.IIP3.V / 2);
    
    dst._channel = FIRFilterConverter!C(model.channel.impulseResponse[0 .. model.channel.taps]);
    dst._rxLNAVGA = PowerControlAmplifierConverter!C((model.thermalNoise.power(model) * (model.INR.dB + model.lna.NF.dB).dB.gain^^2).sqrt.V);
    dst._rxLNARapp = RappModelConverter!C(model.lna.GAIN, model.lna.smoothFactor, (model.lna.IIP3.dBm - 36).dB.gain);
    dst._rxIQMixer = IQImbalanceConverter!C(0.dB, model.rxIQMixer.imbCoef);
    dst._rxQZVGA = PowerControlAmplifierConverter!C((30 - model.ofdm.PAPR.dB + 4.76).dBm);
    dst._rxQZ = SimpleQuantizerConverter!C(model.quantizer.numOfBits);

    return dst;
}
