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

import dffdd.blockdiagram.adder;
import dffdd.blockdiagram.amplifier;
import dffdd.blockdiagram.utils;
import dffdd.utils.unit;

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
            assert(txBaseband.front == txBasebandSWP.front);
            assert(paDirect.front == paDirectSWP.front);

            assert(receivedSI.front == receivedSISWP.front);
            assert(receivedSI.front == received.front);

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


SimulatedSignals makeSimulatedSignals(Model model, string resultDir = null)
{
    immutable doOutput = resultDir !is null;

    bool* swIsAGCTraining = new bool(true),
          swIsNotAGCTraining = new bool(false);

    immutable bool* alwaysFalsePointer = new bool(false);
    immutable bool* alwaysTruePointer = new bool(true);

    SimulatedSignals dst = new SimulatedSignals;
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

    return dst;
}
