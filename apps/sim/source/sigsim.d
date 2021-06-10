module sigsim;

import std.algorithm;
import std.complex;
import std.file;
import std.format;
import std.math;
import std.meta;
import std.path;
import std.random;
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
import dffdd.dpd.polynomial;
import dffdd.utils.unit;
import dffdd.utils.json;

import models;
import simmain;
import snippet;


alias Signal = ForwardRange!(Complex!double);


Signal asSignal(R)(R r)
if(isForwardRange!R)
{
    return r.toWrappedRange;
}


final class SimulatedSignals
{
    alias C = Complex!double;


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

        // まずはDPDを無視して，VGAを学習させる
        while(!this.isConvergedAllVGA([&_txDPD, &_detxDPD])){
            this.popFrontN(_model.ofdm.numOfSamplesOf1Symbol);
        }

        this.popFrontN(_model.ofdm.numOfSamplesOf1Symbol * 4);

        auto paDirect = new C[_model.ofdm.numOfSamplesOf1Symbol * _model.dpd.numOfTrainingSymbols];
        auto txBaseband = new C[_model.ofdm.numOfSamplesOf1Symbol * _model.dpd.numOfTrainingSymbols];
        // 次にDPDの学習をする
        if(! _txDPD.isNull) {
            this.fillBuffer!(["txBaseband", "paDirect"])(txBaseband, paDirect);
            _txDPD.get().estimate(txBaseband, paDirect);
        }

        if(! _detxDPD.isNull) {
            this.fillBuffer!(["desiredBaseband", "desiredPADirect"])(txBaseband, paDirect);
            _detxDPD.get().estimate(txBaseband, paDirect);
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
    do{
        // static foreach(m; aliasSeqOf!ms)
        //     static assert(canFind(["txBaseband", "desiredBaseband", "paDirect", "desiredPADirect", "received"], m), m ~ " is illegal.");

        immutable size_t len = buffers[0].length;
        if(len == 0) return;

        foreach(ref buf; _tempbuf)
            buf.length = max(buf.length, len);

        C[] ds = _tempbuf[0][0 .. len],
            xs = _tempbuf[1][0 .. len],
            ns = _tempbuf[2][0 .. len];

        foreach(i; 0 .. len) {
            ds[i] = desiredBaseband.front;
            xs[i] = txBaseband.front;
            ns[i] = noise.front;

            desiredBaseband.popFront();
            txBaseband.popFront();
            noise.popFront();
        }

        C[] dpdds = _tempbuf[20][0 .. len],
            dpdxs = _tempbuf[21][0 .. len];
        
        if(_model.useSTXDPD) _txDPD.get().apply(xs, dpdxs);
        else                 dpdxs[] = xs[];

        if(_model.useDTXDPD) _detxDPD.get().apply(ds, dpdds);
        else                 dpdds[] = ds[];

        C[] txiqs = _tempbuf[3][0 .. len],
            txvgas = _tempbuf[4][0 .. len],
            txpas = _tempbuf[5][0 .. len],
            detxiqs = _tempbuf[6][0 .. len],
            detxvgas = _tempbuf[7][0 .. len],
            detxpas = _tempbuf[8][0 .. len];

        C[] txpn = _tempbuf[17][0 .. len],
            depn = _tempbuf[18][0 .. len],
            onesbuf = _tempbuf[19][0 .. len];

        onesbuf[] = C(1);
        txpn[] = C(1);
        depn[] = C(1);

        // 自端末の送信機のIQインバランス
        if(_model.useSTXIQ) _txIQMixer.get()(xs, txiqs);
        else                txiqs[] = xs[];

        // 相手端末の送信機のIQインバランス
        if(_model.useDTXIQ) _detxIQMixer.get()(ds, detxiqs);
        else                detxiqs[] = ds[];

        // 自端末と相手端末の位相雑音を生成
        if(_model.useSTXPN) _pnGen1.get()(onesbuf, txpn);
        if(_model.useDTXPN) _pnGen2.get()(onesbuf, depn);

        // 両端末の送信機で位相雑音をかける
        foreach(i; 0 .. len) {
            txiqs[i] *= txpn[i];
            detxiqs[i] *= depn[i];
        }

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

        // 自端末の位相雑音で割る
        foreach(i; 0 .. len) {
            rxlnas[i] /= txpn[i];
        }

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
            else static if(m == "desiredPADirect")
                detxpas.copy(buffers[i]);
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

        if(!this._pnGen1.isNull)        dst._pnGen1 = this._pnGen1.get.dup;
        if(!this._pnGen2.isNull)        dst._pnGen2 = this._pnGen2.get.dup;

        if(!this._txIQMixer.isNull)     dst._txIQMixer = this._txIQMixer.get.dup;
        if(!this._txPAVGA.isNull)       dst._txPAVGA = this._txPAVGA.get.dup;
        if(!this._txPANonlin.isNull)    dst._txPANonlin = cast(IAmplifier!C)this._txPANonlin.get.dup;

        if(!this._detxIQMixer.isNull)     dst._detxIQMixer = this._detxIQMixer.get.dup;
        if(!this._detxPAVGA.isNull)       dst._detxPAVGA = this._detxPAVGA.get.dup;
        if(!this._detxPANonlin.isNull)    dst._detxPANonlin = cast(IAmplifier!C)this._detxPANonlin.get.dup;
        
        if(!this._channelSI.isNull)         dst._channelSI = this._channelSI.get.dup;
        if(!this._channelDesired.isNull)    dst._channelDesired = this._channelDesired.get.dup;

        if(!this._rxDESVGA.isNull)      dst._rxDESVGA = this._rxDESVGA.get.dup;
        if(!this._rxLNAVGA.isNull)      dst._rxLNAVGA = this._rxLNAVGA.get.dup;
        if(!this._rxLNANonlin.isNull)   dst._rxLNANonlin = cast(IAmplifier!C)this._rxLNANonlin.get.dup;
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


    bool isConvergedAllVGA(const(void*)[] ignoreList = null) @property
    {
        bool _isconverged(X)(X* v)
        {
            if(ignoreList.canFind(cast(void*)v))
                return true;

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
    C[][22] _tempbuf;

    bool _nowTrainingMode;
    bool* _useSWPOFDM;

    C _desiredCoef = C(1);
    C _selfInterferenceCoef = C(1);
    C _noiseCoef = C(1);

    Nullable!(FROPhaseNoiseGenerator!(C, Xorshift)) _pnGen1;
    Nullable!(FROPhaseNoiseGenerator!(C, Xorshift)) _pnGen2;

    Nullable!(PolynomialPredistorter!C) _txDPD;
    Nullable!(PolynomialPredistorter!C) _detxDPD;

    Nullable!(IQImbalanceConverter!C) _txIQMixer;
    Nullable!(PowerControlAmplifierConverter!C) _txPAVGA;
    Nullable!(IAmplifier!C) _txPANonlin;

    Nullable!(IQImbalanceConverter!C) _detxIQMixer;
    Nullable!(PowerControlAmplifierConverter!C) _detxPAVGA;
    Nullable!(IAmplifier!C) _detxPANonlin;

    Nullable!(FIRFilterConverter!C) _channelSI;
    Nullable!(FIRFilterConverter!C) _channelDesired;

    Nullable!(PowerControlAmplifierConverter!C) _rxDESVGA;
    Nullable!(PowerControlAmplifierConverter!C) _rxLNAVGA;
    Nullable!(IAmplifier!C) _rxLNANonlin;
    Nullable!(IQImbalanceConverter!C) _rxIQMixer;
    Nullable!(PowerControlAmplifierConverter!C) _rxQZVGA;
    Nullable!(SimpleQuantizerConverter!C) _rxQZ;
}


SimulatedSignals makeSimulatedSignals(Model model, string resultDir = null)
{
    alias C = Complex!double;

    immutable doOutput = resultDir !is null;
    immutable bool* alwaysFalsePointer = new bool(false);
    immutable bool* alwaysTruePointer = new bool(true);

    SimulatedSignals dst = new SimulatedSignals();
    dst._model = model;
    dst._useSWPOFDM = new bool(false);

    dst.desiredBaseband = desiredRandomBits(model)
                        .connectToModulator(modOFDM(model), alwaysFalsePointer, model)
                        .map!(a => cast(C)(a*1.0))
                        .connectTo!(PowerControlAmplifierConverter!C)(30.dBm, 1e-2)
                        .asSignal;

    dst.txBaseband = siRandomBits(model)
                        .connectToModulator(modOFDM(model), dst._useSWPOFDM, model)
                        .map!(a => cast(C)(a*1.0))
                        .connectTo!(PowerControlAmplifierConverter!C)(30.dBm, 1e-2)
                        .asSignal;

    dst.noise = thermalNoise(model).connectTo!VGA(model.lna.NF).map!(a => cast(C)(a*1.0)).toWrappedRange;

    if(model.useSTXPN)
        dst._pnGen1 = FROPhaseNoiseGenerator!(C, Xorshift)(Xorshift( (model.rndSeed + hashOf("STXPN")) & uint.max ), model.samplingFreq, model.phaseNoise.betaBWHz);

    if(model.useDTXPN)
        dst._pnGen2 = FROPhaseNoiseGenerator!(C, Xorshift)(Xorshift( (model.rndSeed + hashOf("DTXPN")) & uint.max ), model.samplingFreq, model.phaseNoise.betaBWHz);

    if(model.useSTXIQ) {
        immutable normgain = 1/sqrt(1 + model.txIQMixer.imbCoef.sqAbs);
        dst._txIQMixer = IQImbalanceConverter!C(Gain.fromVoltageGain(normgain), model.txIQMixer.imbCoef);
    }

    if(model.useDTXIQ) {
        immutable normgain = 1/sqrt(1 + model.txIQMixer.imbCoef.sqAbs);
        dst._detxIQMixer = IQImbalanceConverter!C(Gain.fromVoltageGain(normgain), model.txIQMixer.imbCoef);
    }

    if(model.useSTXPA){
        dst._txPAVGA = PowerControlAmplifierConverter!C(model.pa.TX_POWER / model.pa.GAIN, 1e-2);

        if(model.pa.modelName == "Rapp")
            dst._txPANonlin = RappModelConverter!C(model.pa.smoothFactor, model.pa.GAIN, model.pa.Vsat).toAmplifierObject!C;
        else if(model.pa.modelName == "Saleh")
            dst._txPANonlin = SalehModelConverter!C(model.pa.GAIN, model.pa.Vsat, PI/6).toAmplifierObject!C;
        else if(model.pa.modelName == "Saleh_noPM")
            dst._txPANonlin = SalehModelConverter!C(model.pa.GAIN, model.pa.Vsat, 0).toAmplifierObject!C;
        else
            assert(0, "Invalid model name of PA");
    }else{
        dst._txPAVGA = PowerControlAmplifierConverter!C(model.pa.TX_POWER, 1e-2);
    }

    if(model.useDTXPA) {
        dst._detxPAVGA = PowerControlAmplifierConverter!C(model.pa.TX_POWER / model.pa.GAIN, 1e-2);

        if(model.pa.modelName == "Rapp")
            dst._detxPANonlin = RappModelConverter!C(model.pa.smoothFactor, model.pa.GAIN, model.pa.Vsat).toAmplifierObject!C;
        else if(model.pa.modelName == "Saleh")
            dst._detxPANonlin = SalehModelConverter!C(model.pa.GAIN, model.pa.Vsat, PI/6).toAmplifierObject!C;
        else if(model.pa.modelName == "Saleh_noPM")
            dst._detxPANonlin = SalehModelConverter!C(model.pa.GAIN, model.pa.Vsat, 0).toAmplifierObject!C;
        else
            assert(0, "Invalid model name of PA");
    }else{
        dst._detxPAVGA = PowerControlAmplifierConverter!C(model.pa.TX_POWER, 1e-2);
    }

    dst._channelSI = FIRFilterConverter!C(model.channelSI.impulseResponse[0 .. model.channelSI.taps]);
    dst._channelDesired = FIRFilterConverter!C(model.channelDesired.impulseResponse[0 .. model.channelDesired.taps]);
    
    Voltage receivedSIPower = model.thermalNoise.power(model) * model.lna.NF * model.INR;
    if(model.useSRXLN){
        dst._rxLNAVGA = PowerControlAmplifierConverter!C(receivedSIPower, 1e-2);

        if(model.lna.modelName == "Rapp")
            dst._rxLNANonlin = RappModelConverter!C(model.lna.smoothFactor, model.lna.GAIN, model.lna.Vsat).toAmplifierObject!C;
        else if(model.lna.modelName == "Saleh")
            dst._rxLNANonlin = SalehModelConverter!C(model.lna.GAIN, model.lna.Vsat, PI/6).toAmplifierObject!C;
        else if(model.pa.modelName == "Saleh_noPM")
            dst._rxLNANonlin = SalehModelConverter!C(model.lna.GAIN, model.lna.Vsat, 0).toAmplifierObject!C;
        else
            assert(0, "Invalid model name of LNA");
        //  dst._rxLNANonlin = SoftLimitConverter!C(model.lna.GAIN, (model.lna.IIP3 / 36.dBm).asV);
    }else{
        dst._rxLNAVGA = PowerControlAmplifierConverter!C(receivedSIPower, 1e-2);
        // dst._rxLNANonlin = RappModelConverter!C(model.lna.GAIN, model.lna.smoothFactor, real.infinity);
    }

    if(model.useSRXIQ) {
        immutable normgain = 1/sqrt(1 + model.rxIQMixer.imbCoef.sqAbs);
        dst._rxIQMixer = IQImbalanceConverter!C(Gain.fromVoltageGain(normgain), model.rxIQMixer.imbCoef);
    }

    dst._rxQZVGA = PowerControlAmplifierConverter!C((30 - model.ofdm.PAPR.asdB + 4.76).dBm, 1e-2);

    if(model.useSRXQZ)
        dst._rxQZ = SimpleQuantizerConverter!C(model.quantizer.numOfBits);

    // immutable snScaleOFDM = Gain.fromPowerGain(1.0L * model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling / model.ofdm.numOfSubcarrier);
    // immutable noisePowerPerSubcarrier = model.thermalNoise.power(model) * model.lna.NF / snScaleOFDM;

    dst._rxDESVGA = PowerControlAmplifierConverter!C(model.thermalNoise.power(model) * model.lna.NF * model.SNR, 1e-2);

    if(model.useSTXDPD)
        dst._txDPD = PolynomialPredistorter!C(model.dpd.order);

    if(model.useDTXDPD)
        dst._detxDPD = PolynomialPredistorter!C(model.dpd.order);

    return dst;
}
