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
import std.variant;
import std.range;

import dffdd.blockdiagram.adder;
import dffdd.blockdiagram.amplifier;
import dffdd.blockdiagram.utils;
import dffdd.utils.unit;

import models;
import simmain;
import snippet;

import rx : doSubscribe, Disposable;


alias Signal = ForwardRange!(Complex!real);


Signal asSignal(R)(R r)
if(isForwardRange!R)
{
    return r.toWrappedRange;
}


/+
final class SimulatedSignals
{
    // this() {}
    alias Signals = getSymbolsByUDA!(SimulatedSignals, signals);


    Complex!real[] front()
    {
        Complex!real[] dst;
        foreach(i, ref e; Signals) if(e !is null) {
            dst ~= e.front;
        }

        return dst;
    }


    bool empty = false;


    void popFront()
    {
        foreach(i, ref e; Signals) if(e !is null)
            e.popFront();
    }


    SimulatedSignals save()
    {
        SimulatedSignals dst = new SimulatedSignals;
        foreach(i, ref e; Signals) if(e !is null) {
            // remove "this."
            mixin("dst." ~ Signals[i].stringof[5 .. $]) = e.save();
        }

        return dst;
    }


    void fillBuffer(alias ms, C, size_t N)(C[][N] buffers...)
    if(ms.length == N)
    in{
        auto s = buffers[0].length;
        foreach(e; buffers) assert(e.length == s);
    }
    body{
        foreach(i; 0 .. buffers[0].length){
            foreach(j, m; aliasSeqOf!ms)
                mixin(format(`buffers[j][i] = cast(C)this.%s.front;`, m));
            
            this.popFront();
        }
    }


    void trainAGC()
    {
        *_swIsAGCTraining = true;
        *_swIsNotAGCTraining = false;

        foreach(i; 0 .. _model.numOfModelTrainingSymbols * _model.ofdm.numOfSamplesOf1Symbol)
        {
            import dffdd.utils.linalg : approxEqualCpx;

            try{
                if(txBaseband !is null && txBasebandSWP !is null)   assert(txBaseband.front == txBasebandSWP.front);
                if(paDirect !is null && paDirectSWP !is null)       assert(paDirect.front == paDirectSWP.front);

                if(receivedSI !is null && receivedSISWP !is null)   assert(receivedSI.front == receivedSISWP.front);
                if(receivedSI !is null && received !is null)        assert(receivedSI.front == received.front);
            }catch(Throwable o){
                import std.stdio;
                writeln(i);
                if(txBaseband !is null && txBasebandSWP !is null)   writefln("txBaseband - txBasebandSWP: %s",  txBaseband.front - txBasebandSWP.front);
                if(paDirect !is null && paDirectSWP !is null)       writefln("paDirect - paDirectSWP:     %s",  paDirect.front - paDirectSWP.front);

                if(receivedSI !is null && receivedSISWP !is null)   writefln("receivedSI - receivedSISWP: %s", receivedSI.front - receivedSISWP.front);
                if(receivedSI !is null && received !is null)        writefln("receivedSI - received:      %s", receivedSI.front - received.front);
                throw o;
            }

            this.popFront();
        }

        *_swIsAGCTraining = false;
        *_swIsNotAGCTraining = true;

        foreach(i; 0 .. _model.numOfModelTrainingSymbols * _model.ofdm.numOfSamplesOf1Symbol)
            this.popFront();
    }


    void popFrontNSym(size_t n)
    {
        foreach(i; 0 .. n * _model.ofdm.numOfSamplesOf1Symbol)
            this.popFront();
    }


  private:
    Model _model;
    bool* _swIsAGCTraining;
    bool* _swIsNotAGCTraining;

    static struct signals {}

  public @signals:
    Signal desiredBaseband,
           txBaseband,
           txBasebandSWP,
           desiredPADirect,
           paDirect,
           paDirectSWP,
           receivedDesired,
           receivedSI,
           receivedSISWP,
           received,
           noise;
}
+/


class SimulatedBlock(C)
    : IMNSPBlock!(
        TGroup!(
            C,  /* self interference */
            C,  /* desired */
        ),
        TGroup!C, true)
{
    ISPBlock!(C, C, true) sTXBB, dTXBB;
    ISPBlock!(C, C, true) sTXIQ, dTXIQ;
    ISPBlock!(C, C, true) sTXPN, dTXPN;
    ISPBlock!(C, C, true) sTXPA, dTXPA;
    ISPBlock!(C, C, true) sTXCH, dTXCH;
    IMNSPBlock!(TGroup!(C, C), TGroup!C, true) sRXANT;
    ISPBlock!(C, C, true) sRXLN;
    ISPBlock!(C, C, true) sRXIQ;
    ISPBlock!(C, C, true) sRXQZ;
    Disposable[Tuple!(string, string)] disposables;


    this(Model model)
    {
        import std.string : toUpper;

        sTXBB = makeNopBlock!C();
        dTXBB = makeNopBlock!C();
        sRXANT = makeZipAddBlock!(C, 2)();

        foreach(block; AliasSeq!("sTXIQ", "dTXIQ", "sTXPN", "dTXPN", "sTXPA", "dTXPA",
                                 "sTXCH", "dTXCH", "sRXLN", "sRXIQ", "sRXQZ"))
        {
            mixin(q{
                if(model.use%2$s)
                    this.%1$s = make%2$sBlock!C(model);
                else
                    this.%1$s = makeNopBlock!C();
            }.format(block, block.toUpper));
        }

        fullConnect();
    }


    Sinks!(C, C) sinks() @property
    {
        Sinks!(C, C) dst;
        dst[0] = sTXBB.sinks[0];
        dst[1] = dTXBB.sinks[0];

        return dst;
    }


    Sources!C sources() @property
    {
        Sources!C dst;
        dst[0] = sRXQZ.sources[0];
        return dst;
    }


    SimulatedBlock!C dup() @property
    {
        SimulatedBlock!C mod = new SimulatedBlock!C();
        mod.sTXBB = cast(typeof(this.sTXBB))this.sTXBB.dup;
        mod.dTXBB = cast(typeof(this.dTXBB))this.dTXBB.dup;
        mod.sTXIQ = cast(typeof(this.sTXIQ))this.sTXIQ.dup;
        mod.dTXIQ = cast(typeof(this.dTXIQ))this.dTXIQ.dup;
        mod.sTXPN = cast(typeof(this.sTXPN))this.sTXPN.dup;
        mod.dTXPN = cast(typeof(this.dTXPN))this.dTXPN.dup;
        mod.sTXPA = cast(typeof(this.sTXPA))this.sTXPA.dup;
        mod.dTXPA = cast(typeof(this.dTXPA))this.dTXPA.dup;
        mod.sTXCH = cast(typeof(this.sTXCH))this.sTXCH.dup;
        mod.dTXCH = cast(typeof(this.dTXCH))this.dTXCH.dup;
        mod.sRXANT = cast(typeof(this.sRXANT))this.sRXANT.dup;
        mod.sRXLN = cast(typeof(this.sRXLN))this.sRXLN.dup;
        mod.sRXIQ = cast(typeof(this.sRXIQ))this.sRXIQ.dup;
        mod.sRXQZ = cast(typeof(this.sRXQZ))this.sRXQZ.dup;

        mod.fullConnect();
        return mod;
    }


    void fullConnect()
    {
        // dispose all existing connection
        foreach(k, d; disposables)
            d.dispose();

        // reconnect
        disposables[tuple("sTXBB", "sTXIQ")]
            = sTXBB.sources[0].doSubscribe(sTXIQ.sinks[0]);
        disposables[tuple("sTXIQ", "sTXPN")]
            = sTXIQ.sources[0].doSubscribe(sTXPN.sinks[0]);
        disposables[tuple("sTXPN", "sTXPA")]
            = sTXPN.sources[0].doSubscribe(sTXPA.sinks[0]);
        disposables[tuple("sTXPA", "sTXCH")]
            = sTXPA.sources[0].doSubscribe(sTXCH.sinks[0]);

        disposables[tuple("dTXBB", "dTXIQ")]
            = dTXBB.sources[0].doSubscribe(dTXIQ.sinks[0]);
        disposables[tuple("dTXIQ", "dTXPN")]
            = dTXIQ.sources[0].doSubscribe(dTXPN.sinks[0]);
        disposables[tuple("dTXPN", "dTXPA")]
            = dTXPN.sources[0].doSubscribe(dTXPA.sinks[0]);
        disposables[tuple("dTXPA", "dTXCH")]
            = dTXPA.sources[0].doSubscribe(dTXCH.sinks[0]);

        disposables[tuple("sTXCH", "sRXANT")]
            = sTXCH.sources[0].doSubscribe(sRXANT.sinks[0]);
        disposables[tuple("dTXCH", "sRXANT")]
            = dTXCH.sources[0].doSubscribe(sRXANT.sinks[1]);

        disposables[tuple("sRXANT", "sRXLN")]
            = sRXANT.sources[0].doSubscribe(sRXLN.sinks[0]);
        disposables[tuple("sRXLN", "sRXIQ")]
            = sRXLN.sources[0].doSubscribe(sRXIQ.sinks[0]);
        disposables[tuple("sRXIQ", "sRXQZ")]
            = sRXIQ.sources[0].doSubscribe(sRXQZ.sinks[0]);
    }


    bool hasGetter(string key)
    {
        // auto ns = key.findSplit(".");
        // if(!ns[1].empty){
        //     if(auto p = ns[0] in blocks){
        //         return p.hasGetter(ns[2]);
        //     }
        // }

        // return false;
        return false;
    }


    bool hasSetter(string key)
    {
        return false;
    }


    void opIndexAssignImpl(Variant value, string key)
    {
        {
            import core.exception;
            onRangeError();
        }
    }


    void opIndexAssign(T)(T value, string key)
    {
        this.opIndexAssignImpl(Variant(value), key);
    }


    Variant opIndex(string key)
    {
        {
            import core.exception;
            onRangeError();

            Variant null_;
            return null_;
        }
    }


  private:
    this() {}
}


class SimulatedSignals(C, R1, R2)
{
    alias C = Complex!real;


    this(SimulatedBlock!C simblock, R1 siSignal, R2 desired)
    {
        _simblock = simblock;
        _si = siSignal;
        _desired = desired;
    }


    void fillBuffer(alias ms, X, size_t N)(X[][N] buffers...)
    if(ms.length == N)
    in{
        auto s = buffers[0].length;
        foreach(e; buffers) assert(e.length == s);
    }
    body{
        import std.string : toUpper;
        import rx : Disposable, doSubscribe;

        Disposable[ms.length] disp;
        foreach(i, m; aliasSeqOf!ms)
            disp[i] = __traits(getMember, _simblock, m).sources[0].doSubscribe(buffers[i]);

        scope(exit)
            foreach(d; disp) d.dispose();

        foreach(k, ref signal; AliasSeq!(_si, _desired)){
            _inputbuffer.length = buffers[0].length;
            foreach(i; 0 .. _inputbuffer.length){
                assert(!signal.empty);
                _inputbuffer[i] = signal.front;
                signal.popFront();
            }

            // auto b = cast(ISPBlock!(C, C, true))_simblock.blocks[_simblock.txName[k]];
            // b.put(_inputbuffer);
            auto b = [_simblock.sTXBB, _simblock.dTXBB][k];
            b.put(_inputbuffer);
        }
    }



    void popFrontN(size_t n)
    body{
        if(_inputbuffer.length < 64)
            _inputbuffer.length = 64;

        foreach(ref j; 0 .. n){
            foreach(k, signal; AliasSeq!(_si, _desired)){
                foreach(i; 0 .. _inputbuffer.length){
                    assert(!signal.empty);
                    _inputbuffer[i] = signal.front;
                    signal.popFront();
                }

                auto b = [_simblock.sTXBB, _simblock.dTXBB][k];
                b.put(_inputbuffer);
            }

            j += _inputbuffer.length;
        }
    }


  private:
    SimulatedBlock!C _simblock;
    R1 _si;
    R2 _desired;
    C[] _inputbuffer;
}


auto makeSimulatedSignals(C, R1, R2)(SimulatedBlock!C block, R1 si, R2 desired)
{
    return new SimulatedSignals!(C, R1, R2)(block, si, desired);
}


void noLoadRun(C, R1, R2)(SimulatedBlock!C block, R1 range1, R2 range2, size_t n, size_t chunkSize = 0)
{
    auto signals = makeSimulatedSignals(block, range1, range2);
    signals._inputbuffer.length = chunkSize;
    signals.popFrontN(n);
}


// auto makeSimulatedSignal(C, SimModel)()


auto makeNopBlock(C)()
{
    NopConverter!C impl;
    return makeRxBlock(impl);
}


auto makeSimulatedBlock(Model model)
{
    return new SimulatedBlock!(Complex!real)(model);
}

/+
SimulatedBlock!(Complex!real) makeSimulatedBlock(Model model)
{
    import std.conv : to;

    alias C = Complex!real;
    SimulatedBlock!C dst = new SimulatedBlock!C;

    foreach(sd; aliasSeqOf!(["S", "D"]))
    {
        string preBlock;
        foreach(i, blockName; aliasSeqOf!([sd~":TX:IQ", sd~":TX:PN", sd~":TX:PA", sd~":TX:CH"])){
            mixin(q{
                if(model.use%1$s)
                    dst.blocks[blockName] = make%1$sBlock!C(model);
                else
                    dst.blocks[blockName] = makeNopBlock!C();

                if(!preBlock.empty)
                    dst.connectAtoB(preBlock, blockName);
            }.format(blockName.filter!"a != ':'".array.to!string));
        }
    }

    {
        auto rxAntenna = makeZipAddBlock!(C, 2)();
        auto disp1 = (cast(ISPBlock!(C, C, true))dst.blocks["S:TX:CH"]).observableObject.doSubscribe(rxAntenna.sinks[0]);
        auto disp2 = (cast(ISPBlock!(C, C, true))dst.blocks["D:TX:CH"]).observableObject.doSubscribe(rxAntenna.sinks[1]);
    
        dst.disposables[tuple("S:TX:CH", "S:RX:ANT")] = disp1;
        dst.disposables[tuple("D:TX:CH", "S:RX:ANT")] = disp2;
    }

    return dst;

    /*
    immutable doOutput = resultDir !is null;

    bool* swIsAGCTraining = new bool(true),
          swIsNotAGCTraining = new bool(false);

    immutable bool* alwaysFalsePointer = new bool(false);
    immutable bool* alwaysTruePointer = new bool(true);

    SimulatedBlock!C dst = new SimulatedBlock!C;
    dst._model = model;
    dst._swIsAGCTraining = swIsAGCTraining;
    dst._swIsNotAGCTraining = swIsNotAGCTraining;

    dst.desiredBaseband = desiredRandomBits(model)
                        .connectToModulator(modOFDM(model), alwaysFalsePointer, model)
                        .map!"a*1.0L".asSignal;

    Signal makeDesired(Signal desiredBaseband)
    {
        // 所望信号
        Signal desired = desiredBaseband;

        if(doOutput){
            desired = desired.psdSaveTo("psd_desired_afterMD.csv", resultDir,
                                        model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        }

        if(model.useDTXIQ){
            desired = desired.connectToTXIQMixer(model);
            if(doOutput){
                desired = desired.psdSaveTo("psd_desired_afterTXIQ.csv", resultDir,
                                            model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
            }
        }

        // if(model.useDTXPN){
        //     desired = desired.connectToTXIQPhaseNoise(model).toWrappedRange;
        //     if(doOutput){
        //         desired = desired.psdSaveTo("psd_desired_afterTXIQPhaseNoise.csv", resultDir,
        //                                     model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        //     }
        // }

        if(model.useDTXPA){
            desired = desired.connectToPowerAmplifier(model);
            if(doOutput){
                desired = desired.psdSaveTo("psd_desired_afterPA.csv", resultDir,
                                            model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
            }
        }

        return desired;
    }

    auto desired = makeDesired(dst.desiredBaseband.save);
    dst.desiredPADirect = desired.save;

    desired = desired.connectTo!PowerControlAmplifier((model.thermalNoise.power(model) * (model.SNR.dB + model.lna.NF.dB).dB.gain^^2).sqrt.V)
                    .toWrappedRange;
    dst.receivedDesired = desired.save;

    if(doOutput){
        desired = desired.psdSaveTo("psd_desired_afterVGA.csv", resultDir,
                                    model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
    }


    Signal makeTXReplica(const(bool)* sw, bool saveToFile)
    {
        Signal txbaseband = siRandomBits(model)
                            .connectToModulator(modOFDM(model), sw, model)
                            .map!"a*1.0L".asSignal;

        if(saveToFile){
            txbaseband = txbaseband.psdSaveTo("psd_SI_afterMD.csv", resultDir,
                                              model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        }

        return txbaseband;
    }

    dst.txBaseband = makeTXReplica(alwaysFalsePointer, false);
    dst.txBasebandSWP = makeTXReplica(swIsNotAGCTraining, false);
    auto selfInterference = makeTXReplica(alwaysFalsePointer, doOutput);


    Signal makePADirect(Signal r, bool saveToFile)
    {
        if(model.useSTXIQ){
            r = r.connectToTXIQMixer(model);
            if(saveToFile) r = r.psdSaveTo("psd_SI_afterTXIQ.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        }

        // if(model.useSTXPN){
            // r = r.connectToTXIQPhaseNoise(model).toWrappedRange;
            // if(saveToFile) r = r.psdSaveTo("psd_SI_afterTXIQPhaseNoise.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        // }

        if(model.useSTXPA){
            r = r.connectToPowerAmplifier(model);
            if(saveToFile) r = r.psdSaveTo("psd_SI_afterPA.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        }
        
        if(model.useSTXIQ2){
            r = r.connectToTXIQMixer(model);
            if(saveToFile) r = r.psdSaveTo("psd_SI_afterTXIQ2.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        }

        return r;
    }

    dst.paDirect = makePADirect(dst.txBaseband.save, false);
    dst.paDirectSWP = makePADirect(dst.txBasebandSWP.save, false);
    selfInterference = makePADirect(selfInterference, doOutput);


    Signal makeReceived(Signal r2, Signal d, Signal n, bool saveToFile, bool forNoise = false)
    {
        auto r = n;

        if(model.withSI && r2 !is null){
            r2 = r2.connectToMultiPathChannel(model)
                .connectTo!PowerControlAmplifier((model.thermalNoise.power(model) * (model.INR.dB + model.lna.NF.dB).dB.gain^^2).sqrt.V)
                .toWrappedRange;

            if(saveToFile) r2 = r2.psdSaveTo("psd_SI_afterVGA.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);

            if(forNoise)
                r = r.add(r2.connectToSwitch(swIsAGCTraining)).toWrappedRange;
            else
                r = r.add(r2).toWrappedRange;
        }

        if(d !is null){
            r = r.add(d.connectToSwitch(swIsNotAGCTraining)).toWrappedRange;
        }

        // r = r.add(thermalNoise(model).psdSaveTo("psd_thermal_noise.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model))
        //         .toWrappedRange;
        if(saveToFile) r = r.psdSaveTo("psd_rcv_afterAWGN.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);

        if(model.useSRXLN){
            r = r.connectToLNA(model).psdSaveTo("psd_rcv_afterLNA.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        }

        if(model.useSRXIQ){
            r = r.connectToRXIQMixer(model);
            
            if(saveToFile)
                r = r.psdSaveTo("psd_rcv_afterRXIQ.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        }

        if(model.useSRXQZ)
            r = r.connectToQuantizer(model);

        return r;
    }

    Signal noise = thermalNoise(model).toWrappedRange;
    // if(doOutput)
    //     dst.noise = noise.save.

    dst.receivedSI = makeReceived(dst.paDirect.save, null, noise.save, false);
    dst.receivedSISWP = makeReceived(dst.paDirectSWP.save, null, noise.save, false);
    dst.received = makeReceived(selfInterference, desired.save, noise.save, doOutput);

    dst.noise = makeReceived(dst.paDirect.save, null, noise.save, false, true);

    if(doOutput)
        dst.noise = dst.noise.psdSaveTo("psd_thermal_noise.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);

    if(!model.useDesiredBaseband) dst.desiredBaseband = null;
    if(!model.useTxBaseband) dst.txBaseband = null;
    if(!model.useTxBasebandSWP) dst.txBasebandSWP = null;
    if(!model.useDesiredPADirect) dst.desiredPADirect = null;
    if(!model.usePADirect) dst.paDirect = null;
    if(!model.usePADirectSWP) dst.paDirectSWP  = null;
    if(!model.useReceivedDesired) dst.receivedDesired = null;
    if(!model.useReceivedSI) dst.receivedSI = null;
    if(!model.useReceivedSISWP) dst.receivedSISWP = null;
    if(!model.useReceived) dst.received = null;
    if(!model.useNoise) dst.noise = null;

    return dst;
    */
}
+/