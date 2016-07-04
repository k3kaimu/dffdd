module app;

import core.thread;

import std.random;
import std.algorithm;
import std.range;
import std.conv;
import std.complex;
import std.stdio;
//import msgpackrpc;
import std.datetime;
import std.math;
import std.functional;
import std.mathspecial;
import std.format;
import std.parallelism;
import std.path;
import std.file : mkdirRecurse;
import std.meta;

import dffdd.blockdiagram.utils;
import dffdd.utils.fft;
import dffdd.blockdiagram.amplifier;
import dffdd.utils.unit;
import dffdd.blockdiagram.noise;
import dffdd.blockdiagram.decimator;
import dffdd.blockdiagram.filter;
import dffdd.blockdiagram.adder;
import dffdd.blockdiagram.mod.ofdm;
import dffdd.blockdiagram.iqmixer;
import dffdd.blockdiagram.quantizer;
//import dffdd.utils.msgpackrpc;


import std.stdio;
import std.algorithm;
import std.typecons;
import std.range;
import std.random;

//import carbon.range;
import dranges.range;
import dranges.algorithm;

import dffdd.itpp;
import dffdd.mod.primitives;
import dffdd.mod.qam;
import dffdd.mod.ofdm;

import dffdd.dsp.statistics;


import constant;
import models;
import snippet;


real qfunc(real x) { return 0.5 * erfc(x / SQRT2); }


//alias Constant = ConstantList.ConstExample1;
//mixin Model!Constant;

//enum ptrdiff_t blockSize = 1024;
//enum ptrdiff_t totalBlocks = 400;


auto psdSaveTo(R)(R r, string filename, string resultDir, size_t dropSize, Model model)
{
    return r.tee(makeSpectrumAnalyzer!cfloat(filename, resultDir, dropSize, model)).toWrappedRange;
}


auto makeSpectrumAnalyzer(C)(string filename, string resultDir, size_t dropSize, Model model)
{
    return makeInstrument!C(delegate void(FiberRange!C r){
        r.drop(dropSize).writePSD(File(buildPath(resultDir, filename), "w"), model.samplingFreq, 1024);
    });
}


void mainImpl(Model model, string resultDir)
{
    immutable ofdmModSignalPower = (){
        auto _modOFDMTest = modOFDM(model);
        return randomBits(1, model).connectToModulator(_modOFDMTest, model).measurePower(1024*1024);
    }();

    mkdirRecurse(resultDir);

    bool* switchDS = new bool(false),
          switchSI = new bool(true);

    auto noise = thermalNoise(model);
    auto received = desiredRandomBits(model)
                    .connectToModulator(modOFDM(model), model)
                        .psdSaveTo("psd_desired_afterMD.csv", resultDir, model.numOfModelTrainingSample, model)
                    .connectToTXIQMixer(model)
                        .psdSaveTo("psd_desired_afterTXIQ.csv", resultDir, model.numOfModelTrainingSample, model)
                    .connectToPowerAmplifier(model)
                        .psdSaveTo("psd_desired_afterPA.csv", resultDir, model.numOfModelTrainingSample, model)
                    .connectTo!PowerControlAmplifier((model.thermalNoise.power(model) * model.SNR.dB.gain^^2).sqrt.V)
                        .psdSaveTo("psd_desired_afterVGA.csv", resultDir, model.numOfModelTrainingSample, model)
                    .connectToSwitch(switchDS)
                    .add(
                        // Self-Interference
                        siRandomBits(model)
                        .connectToModulator(modOFDM(model), model)
                            .psdSaveTo("psd_SI_afterMD.csv", resultDir, model.numOfModelTrainingSample, model)
                        .connectToTXIQMixer(model)
                            .psdSaveTo("psd_SI_afterTXIQ.csv", resultDir, model.numOfModelTrainingSample, model)
                        .connectToPowerAmplifier(model)
                            .psdSaveTo("psd_SI_afterPA.csv", resultDir, model.numOfModelTrainingSample, model)
                        .connectTo!PowerControlAmplifier((model.thermalNoise.power(model) * model.INR.dB.gain^^2).sqrt.V)
                            .psdSaveTo("psd_SI_afterVGA.csv", resultDir, model.numOfModelTrainingSample, model)
                        .connectToSwitch(switchSI)
                    )
                    .connectToAWGN(model)
                        .psdSaveTo("psd_rcv_afterAWGN.csv", resultDir, model.numOfModelTrainingSample, model)
                    .connectToLNA(model)
                        .psdSaveTo("psd_rcv_afterLNA.csv", resultDir, model.numOfModelTrainingSample, model)
                    .connectToRXIQMixer(model)
                        .psdSaveTo("psd_rcv_afterRXIQ.csv", resultDir, model.numOfModelTrainingSample, model)
                    .connectToQuantizer(model)
                    ;

    auto txReplica = siRandomBits(model).connectToModulator(modOFDM(model), model);

    // モデルの安定化
    *switchDS = false;
    *switchSI = true;
    received.popFrontN(model.numOfModelTrainingSample);
    txReplica.popFrontN(model.numOfModelTrainingSample);

    auto recvs = new Complex!float[model.blockSize],
         refrs = new Complex!float[model.blockSize],
         outps = new Complex!float[model.blockSize];

    auto fftObj = new Fft(model.blockSize);
    auto filter = makeCascadeHammersteinFilter(modOFDM(model), model);

    {
        File powerFile = File(buildPath(resultDir, "errorout_long.csv"), "w");

        foreach(blockIdx; 0 .. model.numOfFilterTrainingSample / model.blockSize)
        {
            foreach(i; 0 .. model.blockSize){
                recvs[i] = (a => complex(a.re, a.im))(received.front);
                refrs[i] = (a => complex(a.re, a.im))(txReplica.front);

                received.popFront();
                txReplica.popFront();
            }

            filter.apply!true(refrs, recvs, outps);

            if(blockIdx == 0)
            {
                File outFile = File(buildPath(resultDir, "errorout_start.csv"), "w");
                foreach(i; 0 .. model.blockSize){
                    if(i > 20000) break;
                    auto pr = recvs[i].re^^2 + recvs[i].im^^2,
                         po = outps[i].re^^2 + outps[i].im^^2,
                         c = -10*log10(po / pr);

                    outFile.writefln("%s,%s,%s,%s,", i, pr, po, c);
                }
            }

            // fft
            {
                auto outputSpec = fftObj.fftWithSwap(outps);
                auto recvSpec = fftObj.fftWithSwap(recvs);

                real sumR = 0, sumO = 0;
                foreach(i; 0 .. model.blockSize){
                    immutable po = outps[i].re^^2 + outps[i].im^^2,
                              pr = recvs[i].re^^2 + recvs[i].im^^2;

                    immutable freq = i * model.samplingFreq / model.blockSize - (model.samplingFreq/2);
                    if(abs(freq)/model.samplingFreq < (model.ofdm.numOfSubcarrier*1.0)/(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling)
                    && abs(freq)/model.samplingFreq > 1.0/model.ofdm.numOfFFT/model.ofdm.scaleOfUpSampling){
                        sumR += pr;
                        sumO += po;
                    }
                }

                powerFile.writefln("%s,%s,%s,%s,", model.blockSize*blockIdx, sumR, sumO, 10*log10(sumO / sumR));
            }
        }
    }

    //auto recvPSD = makeInstrument!(Complex!float, r => r.map!"a.re + a.im*1i".writePSD(File(buildPath(resultDir, "psd_beforeSIC.csv"), "w"), Constant.samplingFreq, 1024));
    //auto sicPSD = makeInstrument!(Complex!float, r => r.map!"a.re + a.im*1i".writePSD(File(buildPath(resultDir, "psd_afterSIC.csv"), "w"), Constant.samplingFreq, 1024));
    //auto foutPSD = makeInstrument!(Complex!float, r => r.map!"a.re + a.im*1i".writePSD(File(buildPath(resultDir, "psd_filter_output.csv"), "w"), Constant.samplingFreq, 1024));

    auto recvPSD = makeSpectrumAnalyzer!(Complex!float)("psd_beforeSIC.csv", resultDir, 0, model);
    auto sicPSD = makeSpectrumAnalyzer!(Complex!float)("psd_afterSIC.csv", resultDir, 0, model);
    auto foutPSD = makeSpectrumAnalyzer!(Complex!float)("psd_filter_output.csv", resultDir, 0, model);

    foreach(blockIdxo; 0 .. 64)
    {
        foreach(i; 0 .. model.blockSize){
            recvs[i] = (a => complex(a.re, a.im))(received.front);
            refrs[i] = (a => complex(a.re, a.im))(txReplica.front);

            received.popFront();
            txReplica.popFront();
        }

        filter.apply!true(refrs, recvs, outps);

        foreach(i; 0 .. model.blockSize)
            refrs[i] = recvs[i] - outps[i];

        .put(recvPSD, recvs);
        .put(sicPSD, outps);
        .put(foutPSD, refrs);
    }


    // BER count
    *switchDS = true;
    real berResult = -1;
    auto berCounter = makeInstrument(delegate void (FiberRange!cfloat r){
        auto bits= r
        .connectTo!PowerControlAmplifier(ofdmModSignalPower.sqrt.V)
        .connectToDemodulator(modOFDM(model), model);

        auto refBits = desiredRandomBits(model);
        berResult = measureBER(bits, refBits, model.berCounter.totalBits);
    });

    while(berResult == -1){
        foreach(i; 0 .. model.blockSize){
            recvs[i] = (a => complex(a.re, a.im))(received.front);
            refrs[i] = (a => complex(a.re, a.im))(txReplica.front);

            received.popFront();
            txReplica.popFront();
        }

        if(model.withSIC)
            filter.apply!false(refrs, recvs, outps);
        else
            outps[] = recvs[];

        foreach(e; outps){
            berCounter.put(e.re + e.im * 1i);
        }
    }

    File file = File(buildPath(resultDir, "ber.csv"), "w");
    file.writeln(berResult);
}


void main()
{
    Model[] models;
    string[] dirs;

    // no sic
    foreach(p; 0 .. 31){
        Model model;
        model.SNR = p;
        model.INR = 30;
        model.withSIC = false;
        //mainImpl(model, format("no_sic_snr%s_inr%s", model.SNR, model.INR));
        models ~= model;
        dirs ~= "no_sic_snr%s_inr%s";
    }

    // with sic
    foreach(p; 0 .. 31){
        Model model;
        model.SNR = p;
        model.INR = 30;
        model.withSIC = true;
        //mainImpl(model, format("with_sic_snr%s_inr%s", model.SNR, model.INR));
        models ~= model;
        dirs ~= "with_sic_snr%s_inr%s";
    }

    foreach(i; iota(models.length).parallel()){
        mainImpl(models[i], dirs[i]);
    }
}


__EOF__


void main()
{
    auto _modOFDMTest = modOFDM(4);
    immutable ofdmModSignalPower = randomBits(1).connectToModulator(_modOFDMTest).measurePower(1024*1024);

    foreach(pSI; iota(20, 100, 10).array.parallel)
    //foreach(p; iota(10, 21))
    {
        auto resultDir = "snr_%s".format(p);
        mkdirRecurse(resultDir);

        //StopWatch sw;
        //sw.start();
        ulong[40] cntList;
        //immutable p = 23;
        //immutable totalBits = p > 20 ? Constant.BERCounter.totalBits : Constant.BERCounter.totalBits / 10;
        immutable totalBits = Constant.BERCounter.totalBits;

        //writeln(ofdmModSignalPower);

        //foreach(pIdx; iota(40).parallel())
        //{
            auto ofdmDesired = modOFDM(4);
            auto ofdmSI = modOFDM(4);
            //auto ofdmX32 = modOFDM(32);
            //Random r1, r2;
            //auto desiredBits = desiredRandomBits();

            //auto _foo_ = PowerControlAmplifier.makeBlock(thermalNoise(), 1.V);
            //real* switchSIGain = new real,
                  //switchDSGain = new real;
            //*switchSIGain = 1;
            //*switchDSGain = 1;
            bool* switchDS = new bool(false),
                  switchSI = new bool(true);

            auto noise = thermalNoise();
            auto received = desiredRandomBits().connectToModulator(ofdmDesired)
                                .binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(model.trainingSamplesOfModel).writePSD(File(buildPath(resDir, "psd_DS_afterMD.csv"), "w"), Constant.samplingFreq, 1024))).toWrappedRange)(resultDir)
                            .connectToTXIQMixer()
                                .binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(model.trainingSamplesOfModel).writePSD(File(buildPath(resDir, "psd_DS_afterIQ.csv"), "w"), Constant.samplingFreq, 1024))).toWrappedRange)(resultDir)
                            .connectToPowerAmplifier()
                                .binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(model.trainingSamplesOfModel).writePSD(File(buildPath(resDir, "psd_DS_afterPA.csv"), "w"), Constant.samplingFreq, 1024))).toWrappedRange)(resultDir)
                            .connectTo!PowerControlAmplifier((Constant.ThermalNoise.power * pD.dB.gain^^2).sqrt.V)
                                .binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(model.trainingSamplesOfModel).writePSD(File(buildPath(resDir, "psd_DS_afterVGA.csv"), "w"), Constant.samplingFreq, 1024))).toWrappedRange)(resultDir)
                            .connectToSwitch(switchDS)
                            .add(
                                // Self-Interference
                                siRandomBits().connectToModulator(ofdmSI)
                                    .binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(model.trainingSamplesOfModel).writePSD(File(buildPath(resDir, "psd_SI_afterMD.csv"), "w"), Constant.samplingFreq, 1024))).toWrappedRange)(resultDir)
                                .connectToTXIQMixer()
                                    .binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(model.trainingSamplesOfModel).writePSD(File(buildPath(resDir, "psd_SI_afterMD.csv"), "w"), Constant.samplingFreq, 1024))).toWrappedRange)(resultDir)
                                .connectToPowerAmplifier()
                                .connectTo!PowerControlAmplifier((Constant.ThermalNoise.power * pSI.dB.gain^^2).sqrt.V)
                                .connectToSwitch(switchSI)
                            )
                            .connectToAWGN()
                            .connectToLNA()
                            .connectToRXIQMixer()
                            .connectToQuantizer()
                            //.connectToTXChain()
                                //.binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(100_000).writePSD(File(buildPath(resDir, "psd_afterMD.csv"), "w"), Constant.samplingFreq, 1024))).toWrappedRange)(resultDir)
                            //.tmap!"a * (*b)"(switchSIGain)
                                //.binaryFun!((r, idx) => idx == 0 ? (r.tee(makeInstrument!(cfloat, r => r.writePSD(File("psd_afterMD_%sdB_1.csv", "w"), Constant.samplingFreq, 1024))).toWrappedRange) : r.toWrappedRange)(pIdx)
                            //.connectToUpSampler()
                                //.binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(100_000).writePSD(File(buildPath(resDir, "psd_afterUS.csv"), "w"), Constant.samplingFreq * Constant.UpDownSampler.scaleOfUpSampling, 1024 * Constant.UpDownSampler.scaleOfUpSampling))).toWrappedRange)(resultDir)
                                //.binaryFun!((r, idx) => idx == 0 ? (r.tee(makeInstrument!(cfloat, r => r.writePSD(File("psd_afterUS_%sdB_1.csv", "w"), Constant.samplingFreq * Constant.UpDownSampler.scaleOfUpSampling, 1024 * Constant.UpDownSampler.scaleOfUpSampling))).toWrappedRange) : r.toWrappedRange)(pIdx)
                            //.connectToTXIQMixer()
                            //.connectToPowerAmplifier()
                            //.connectTo!VGA((-(10*log10(snr0) - p)).dB)
                            //.connectTo!PowerControlAmplifier((Constant.ThermalNoise.power * p.dB.gain^^2).sqrt.V)
                                //.binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(100_000).writePSD(File(buildPath(resDir, "psd_afterPA.csv"), "w"), Constant.samplingFreq * Constant.UpDownSampler.scaleOfUpSampling, 1024 * Constant.UpDownSampler.scaleOfUpSampling))).toWrappedRange)(resultDir)
                                //.binaryFun!((r, idx) => idx == 0 ? (r.tee(makeInstrument!(cfloat, r => r.writePSD(File("psd_afterPA_%sdB_1.csv", "w"), Constant.samplingFreq * Constant.UpDownSampler.scaleOfUpSampling, 1024 * Constant.UpDownSampler.scaleOfUpSampling))).toWrappedRange) : r.toWrappedRange)(pIdx)
                            //.connectToFlatRayleighFadingChannel()
                            //.connectTo!PowerControlAmplifier((Constant.ThermalNoise.power * p.dB.gain^^2).sqrt.V)
                            //.connectToAWGN()
                                //.binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(100_000).writePSD(File(buildPath(resDir, "psd_afterRX.csv"), "w"), Constant.samplingFreq * Constant.UpDownSampler.scaleOfUpSampling, 1024 * Constant.UpDownSampler.scaleOfUpSampling))).toWrappedRange)(resultDir)
                                //.binaryFun!((r, idx) => idx == 0 ? (r.tee(makeInstrument!(cfloat, r => r.writePSD(File("psd_afterRX_%sdB_1.csv", "w"), Constant.samplingFreq * Constant.UpDownSampler.scaleOfUpSampling, 1024 * Constant.UpDownSampler.scaleOfUpSampling))).toWrappedRange) : r.toWrappedRange)(pIdx)
                            //.connectToLNA()
                            //.connectToRXIQMixer()
                            //.connectToQuantizer()
                            //.connectToDownSampler()
                                //.binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(100_000).writePSD(File(buildPath(resDir, "psd_afterDS.csv"), "w"), Constant.samplingFreq, 1024))).toWrappedRange)(resultDir)
                            //.connectTo!VGA(((10*log10(snr0) - p + 3)).dB)
                            //.connectTo!PowerControlAmplifier(ofdmModSignalPower.sqrt.V)
                                //.binaryFun!((r, resDir) => r.tee(makeInstrument!(cfloat, r => r.drop(100_000).writePSD(File(buildPath(resDir, "psd_afterQZ.csv"), "w"), Constant.samplingFreq, 1024))).toWrappedRange)(resultDir)
                            //.connectToDemodulator(ofdm)
                            ;

            // 定常状態にする
            received.popFrontN(1000*1000);
/+
            auto refBits = siRandomBits();
            //writeln(measureBER(received, refBits, 10_000));
            auto reference = refBits.connectToModulator(ofdmDesired);

            //received.popFrontN(100*1000);
            reference.popFrontN(100*1000);

            //reference.popFrontN(4);

            writeln("Setup is ended");

            Fft fftObj = new Fft(blockSize);

            size_t sumCNT;

            auto recvs = new Complex!float[blockSize],
                 refrs = new Complex!float[blockSize],
                 outps = new Complex!float[blockSize];

            //double[] fftResultRecv = new double[blockSize],
            //         fftResultSIC = new double[blockSize],
            //         fftResultFilterOutput = new double[blockSize];

            //fftResultRecv[] = 0;
            //fftResultSIC[] = 0;
            //fftResultFilterOutput[] = 0;

            File powerFile = File(buildPath(resultDir, "errorout_long.csv"), "w");

            writeln("Some env is setup-ed");

            auto filter = makeCascadeHammersteinFilter(ofdm);
            //auto filter = makeParallelHammersteinFilter(ofdm);

            writeln("Filter is setup-ed");

            foreach(blockIdx; 0 .. totalBlocks)
            {
                foreach(i; 0 .. blockSize){
                    recvs[i] = (a => complex(a.re, a.im))(received.front);
                    refrs[i] = (a => complex(a.re, a.im))(reference.front);

                    received.popFront();
                    reference.popFront();
                }

                filter.apply!true(refrs, recvs, outps);

                // fft
                {
                    auto outputSpec = fftObj.fftWithSwap(outps);
                    auto recvSpec = fftObj.fftWithSwap(recvs);

                    real sumR = 0, sumO = 0;
                    foreach(i; 0 .. blockSize){
                        immutable po = outps[i].re^^2 + outps[i].im^^2,
                                  pr = recvs[i].re^^2 + recvs[i].im^^2;

                        //fftResultRecv[i] += pr;
                        //fftResultSIC[i] += po;
                        //fftResultFilterOutput[i] += (x => x.re^^2 + x.im^^2)(recvs[i] - outps[i]);

                        immutable freq = i * Constant.samplingFreq / blockSize - (Constant.samplingFreq/2);
                        if(abs(freq)/Constant.samplingFreq < (Constant.OFDM.numOfSubcarrier*1.0)/(Constant.OFDM.numOfFFT * Constant.OFDM.scaleOfUpSampling)
                        && abs(freq)/Constant.samplingFreq > 1.0/Constant.OFDM.numOfFFT/Constant.OFDM.scaleOfUpSampling){
                            sumR += pr;
                            sumO += po;
                        }
                    }

                    powerFile.writefln("%s,%s,%s,%s,", blockSize*blockIdx, sumR, sumO, 10*log10(sumO / sumR));
                    ++sumCNT;
                }


                if(blockIdx == 0)
                {
                    File outFile = File(buildPath(resultDir, "errorout_start.csv"), "w");
                    foreach(i; 0 .. blockSize){
                        if(i > 20000) break;
                        auto pr = recvs[i].re^^2 + recvs[i].im^^2,
                             po = outps[i].re^^2 + outps[i].im^^2,
                             c = -10*log10(po / pr);

                        outFile.writefln("%s,%s,%s,%s,", i, pr, po, c);
                    }
                }

                //if(blockIdx == 300*1024/blockSize)
                //{
                //    File outFile = File(buildPath(resultDir, "fftout.csv"), "w");

                //    real sumR = 0, sumO = 0;
                //    foreach(i; 0 .. blockSize)
                //    {
                //        auto pr = fftResultRecv[i],
                //             po = fftResultSIC[i],
                //             before = 10*log10(fftResultRecv[i] / sumCNT),
                //             after = 10*log10(fftResultSIC[i] / sumCNT),
                //             filterOutput = 10*log10(fftResultFilterOutput[i] / sumCNT);

                //        real freq = i * Constant.samplingFreq / blockSize - (Constant.samplingFreq / 2);

                //        if(abs(freq)/Constant.samplingFreq < (Constant.OFDM.numOfSubcarrier*1.0)/(Constant.OFDM.numOfFFT * Constant.OFDM.scaleOfUpSampling)
                //        && abs(freq)/Constant.samplingFreq > 1.0/Constant.OFDM.numOfFFT/Constant.OFDM.scaleOfUpSampling){
                //            sumR += pr;
                //            sumO += po;
                //        }

                //        outFile.writefln("%s,%s,%s,%s,%s,%s,%s", blockIdx, i, freq, filterOutput, before, after, before - after);
                //    }

                //    writefln("%s [dB](%s / %s = %s)", 10*log10(sumO / sumR), sumR, sumO, sumR / sumO);
                //}

                //if(blockIdx % (100 * 1024 / blockSize) == 0){
                //    fftResultRecv[] = 0;
                //    fftResultSIC[] = 0;
                //    fftResultFilterOutput[] = 0;

                //    sumCNT = 0;
                //}
            }

            {
                StopWatch sw;
                sw.start();

                foreach(blockIdx; 0 .. 10000)
                {
                    foreach(i; 0 .. blockSize){
                        recvs[i] = (a => complex(a.re, a.im))(received.front);
                        refrs[i] = (a => complex(a.re, a.im))(reference.front);

                        received.popFront();
                        reference.popFront();
                    }

                    filter.apply!true(refrs, recvs, outps);
                }

                sw.stop();
                writefln("%s[ksps]", 10000 * blockSize * 1.0 / sw.peek.msecs);
            }

            auto recvPSD = makeInstrument!(Complex!float, r => r.map!"a.re + a.im*1i".writePSD(File(buildPath(resultDir, "psd_beforeSIC.csv"), "w"), Constant.samplingFreq, 1024));
            auto sicPSD = makeInstrument!(Complex!float, r => r.map!"a.re + a.im*1i".writePSD(File(buildPath(resultDir, "psd_afterSIC.csv"), "w"), Constant.samplingFreq, 1024));
            auto foutPSD = makeInstrument!(Complex!float, r => r.map!"a.re + a.im*1i".writePSD(File(buildPath(resultDir, "psd_filter_output.csv"), "w"), Constant.samplingFreq, 1024));

            foreach(blockIdxo; 0 .. 64)
            {
                foreach(i; 0 .. blockSize){
                    recvs[i] = (a => complex(a.re, a.im))(received.front);
                    refrs[i] = (a => complex(a.re, a.im))(reference.front);

                    received.popFront();
                    reference.popFront();
                }

                filter.apply!true(refrs, recvs, outps);

                foreach(i; 0 .. blockSize)
                    refrs[i] = recvs[i] - outps[i];

                .put(recvPSD, recvs);
                .put(sicPSD, outps);
                .put(foutPSD, refrs);
            }


            auto noisePSD = makeInstrument!(Complex!float, r => r.map!"a.re + a.im*1i".writePSD(File(buildPath(resultDir, "psd_noise.csv"), "w"), Constant.samplingFreq, 1024));

            *switchSIGain = 0;
            received.popFrontN(blockSize * 4);
            foreach(blockIdx; 0 .. 64)
            {
                foreach(i; 0 .. blockSize){
                    recvs[i] = (a => complex(a.re, a.im))(received.front);
                    received.popFront();
                }
                .put(noisePSD, recvs);
            }

            /*
            ulong cnt;
            foreach(i; 0 .. totalBits / 40){
                auto a = demodBits.front;
                auto b = refBits.front;

                if(a != b) ++cnt;

                //writef("(%s, %s), ", a, b);
                demodBits.popFront();
                refBits.popFront();
            }
            cntList[pIdx] = cnt;
            */
        //}

        //writefln("%s,%e,%f",
        //        p,
        //        (cntList[].sum())*1.0 / totalBits,
        //        totalBits * 1.0 / sw.peek.msecs * 1000
        //);

/*
        auto psd = demodBits.calculatePowerSpectralDensity(Constant.samplingFreq*8, 1024);
        File file = File(format("psd_SNR_%s_2.csv", p), "w");
        foreach(i, e; psd)
            file.writefln("%f,%f,%f", (i*1.0/1024*Constant.samplingFreq*8-(Constant.samplingFreq*8/2))/1e6, e, e == 0 ? -400 : 30+10*log10(e));
*/
    +/
    }

/*
    auto fftObj = new Fft(1024);
    Complex!float[] buf = new Complex!float[1024];
    Complex!float[] fftRes = new Complex!float[1024];
    buf[] = complex!float(0, 0);
    foreach(i; 0 .. 1024){
        real x2 = (cast(real)i) - 1024/2;
        if(x2 == 0)
            buf[i] = complex!float(1, 0);
        else
            buf[i] = complex!float(sin(PI*x2/2.0)/(x2*PI/2.0)/2.0, 0);
    }
*/
    //writeln(buf);
    //buf[1] = complex!float(sin(PI * 1)/1/PI/2, 0);
    //buf[3] = complex!float(sin(PI * 3)/3/PI/2, 0);
    //fftRes[] = complex!float(1, 0);
    //fftRes[$/4 .. $-$/4] = complex!float(0, 0);
    //fftObj.inverseFft(fftRes, buf);
    //fftObj.fft(buf, fftRes);
    //writeln(buf);
    //buf[$/2 .. $]

    //writefln(buf);
    //foreach(i, e; fftRes)
        //writefln("%s,%e,%e,%e", i, e == 0 ? -400 : 10*log10(e.re^^2 + e.im^^2), e.arg/PI, 10*log10(buf[i].re^^2));

    //foreach(p; 0 .. 20)
    //{
    //    Random r1;

    //    auto noise = ThermalNoise(SamplingFrequency(samplingFreq*4), 300);
    //    auto signal = r1.map!"cast(ubyte)(a&1)".splitN(ofdmX16.symInputLength).tmap!(reverseArgs!mod, [0])(ofdmX16).joiner
    //                    .connectTo!VGA((-(10*log10(snr0) - p)).dB)
    //                    .zip(noise).map!"cast(cfloat)(a[0]+a[1])";

    //    auto psd = signal.calculatePowerSpectralDensity(samplingFreq*4, 1024);
    //    File file = File(format("psd_SNR_%s.csv", p), "w");
    //    foreach(i, e; psd)
    //        file.writefln("%f,%f,%f", (i*1.0/1024*samplingFreq-(samplingFreq/2))/1e6, e, e == 0 ? -400 : 30+10*log10(e));
    //}

    //writeln(10*log10(signalPower / noisePower));

    //Random r1, r2;


    //auto ofdm = chainedMod(dffdd.mod.qam.QAM(16), dffdd.mod.ofdm.OFDM(64, 16, 48, 1));
    //auto ofdm = setGlobal!"ofdm"(chainedMod(dffdd.mod.qam.QAM(16), dffdd.mod.ofdm.OFDM(64, 16, 48, 1)));

    //writeln(ofdm.symInputLength);
    //writeln(ofdm.symOutputLength);

    //Random r1, r2;

    //auto sig = r1.map!"cast(ubyte)(a&1)"
                //.splitN(ofdm.symInputLength).map!array.tmap!(reverseArgs!mod, [0])(ofdm).joiner()
                //.zip(y.connectTo!VGA(200.dB)).map!"cast(cfloat)(a[0]+a[1])"
                //.splitN(ofdm.symOutputLength).map!array.tmap!(reverseArgs!demod, [0])(ofdm).joiner().take(100);

    //writeln(sig);
    //writeln(r2.map!"cast(ubyte)(a&1)".zip(sig).map!"[a.tupleof]");

    //auto fftRes = calculatePowerSpectralDensity(y, samplingFreq, 1024);
    //size_t i = 0;
    //foreach(e; fftRes)
    //{
    //    writefln("%f,%f,%f", (i*1.0/1024*samplingFreq-(samplingFreq/2))/1e6, e, e == 0 ? -400 : 30+10*log10(e));
    //    ++i;
    //}
}



/+
struct OFDMModulator
{
    static
    auto makeBlock(R)(R r, TCPClient client)
    if(isInputRange!R && is(ElementType!R == bool))
    {
        import std.typecons;

        return RefCounted!(OFDMModulatorImpl!R)(r, client);
    }


    static
    struct OFDMModulatorImpl(R)
    {
        this(R r, TCPClient client)
        {
            _r = r;
            _client = client;
            _empty = false;
            _i = 127;
            _buf.length = 128;
            popFront();
        }


        cfloat front() const @property { return _buf[_i]; }
        bool empty() const @property { return _empty; }

        void popFront()
        {
            ++_i;

            if(_i == _buf.length){
                _i = 0;
                _buf.length = 0;
                generate();
            }
        }

      private:
        R _r;
        TCPClient _client;
        bool _empty;
        size_t _i;
        cfloat[] _buf;


        void generate()
        {
            bool sent;
            size_t cnt;
            while(1){
                auto rep = _client.call!(float[][])("pullOutputs1", 128);
                if(rep.empty && !sent){
                    uint[128] _arr_;
                    uint[] bits = _arr_[];

                    foreach(i; 0 .. 129){
                        if(_r.empty || i == 128){
                            bits = bits[0 .. i];
                            break;
                        }

                        foreach(j; 0 .. 8){
                            bits[i] <<= 1;
                            if(!_r.empty){
                                bits[i] = bits[i] | (_r.front ? 1 : 0);
                                _r.popFront();
                            }
                        }
                    }

                    if(bits.length == 0){
                        assert(_r.empty);
                        ++cnt;
                        if(cnt == 1024){
                            _empty = true;
                            return;
                        }
                        Thread.sleep(dur!"msecs"(10));
                        continue;
                    }
                    _client.notify("pushInputs1", bits);
                }
                else if(rep.empty && sent){
                    break;
                }else{
                    sent = true;
                    foreach(i; 0 .. rep.length)
                        _buf ~= rep[i][0] + rep[i][1]*1i;
                }
            }
        }
    }
}


struct OFDMDemodulator
{
    static
    auto makeBlock(R)(R r, TCPClient client)
    if(isInputRange!R && is(ElementType!R : creal))
    {
        import std.typecons;

        return RefCounted!(OFDMDemodulatorImpl!R)(r, client);
    }


    static
    struct OFDMDemodulatorImpl(R)
    {
        this(R r, TCPClient client)
        {
            _r = r;
            _client = client;
            _empty = false;
            _i = 127;
            _j = 7;
            _buf.length = 128;
            popFront();
        }


        bool front() const @property { return (_buf[_i] & (0b10000000 >> (_j))) > 0; }
        bool empty() const @property { return _empty; }

        void popFront()
        {
            ++_j;
            if(_j == 8){
                ++_i;
                _j = 0;
            }

            if(_i == _buf.length){
                _i = 0;
                _buf.length = 0;
                generate();
            }
        }

      private:
        R _r;
        TCPClient _client;
        bool _empty;
        size_t _i, _j;
        ubyte[] _buf;


        void generate()
        {
            bool sent;
            size_t cnt;
            while(1){
                auto rep = _client.call!(uint[])("pullOutputs2", 128);
                if(rep.empty && !sent){
                    float[2][128] _arr_;
                    float[2][] sigs = _arr_[];

                    foreach(i; 0 .. 129){
                        if(_r.empty || i == 128){
                            sigs = sigs[0 .. i];
                            break;
                        }

                        if(!_r.empty){
                            auto e = _r.front;
                            _r.popFront();

                            sigs[i][0] = e.re;
                            sigs[i][1] = e.im;
                        }
                    }

                    if(sigs.length == 0){
                        assert(_r.empty);
                        ++cnt;
                        if(cnt == 1024){
                            _empty = true;
                            return;
                        }
                        Thread.sleep(dur!"msecs"(10));
                        continue;
                    }
                    _client.notify("pushInputs2", sigs);
                }
                else if(rep.empty && sent){
                    break;
                }else{
                    sent = true;
                    foreach(i; 0 .. rep.length)
                        _buf ~= cast(ubyte)rep[i];
                }
            }
        }
    }
}


enum real samplingFreq = 40e6;
enum size_t upSamplingRatio = 1;
enum real upSampledFreq = upSamplingRatio * samplingFreq;
enum real ncoFreq = 8e6;
enum size_t Res = 1024;


void main()
{
    auto client = new TCPClient(Endpoint(18800, "127.0.0.1"));
    auto y = repeat(1).map!(a => uniform01() < 0.5).cycle.connectTo!(OFDMModulator)(client)
            .connectTo!VGA((10).dB)
            .connectTo!VGA((-30).dB)
            .connectTo!PowerAmplifier(30.dB, 10.dBm);
    //auto y = ThermalNoise(SamplingFrequency(samplingFreq), 300);
    //writeln(y.take(10));

    auto p = calculateAveragePower(y, samplingFreq, 1024);

    writeln(30+10*log10(p));

    foreach(i; 0 .. Res)
        y.popFront();

    auto fftRes = y.calculatePowerSpectralDensity(Res, 32);

    size_t i = 0;
    foreach(e; fftRes)
    {
        writefln("%f,%f,%f", (i*1.0/Res*samplingFreq-(samplingFreq/2))/1e6, e, e == 0 ? -400 : 30+10*log10(e));
        ++i;
    }

    //auto r2 = r.connectTo!(OFDMDemodulator)(client);
    //auto r = OFDMModulator.makeBlock(repeat(1).map!(a => uniform01()).map!"a < 0.5"(), client);
    //writeln(r2);

    /+
    SysTime start = Clock.currTime;
    size_t cnt;

    while(1){
        uint[] demod;
        // mod
        while(demod.empty)
        {
            LLoop:
            while(1){
                client.notify("pushInputs1", [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]);
                //Thread.sleep(dur!"msecs"(100));

                while(1){
                    auto rep = client.call!(float[][])("pullOutputs1", 128);
                    //Thread.sleep(dur!"msecs"(100));
                    if(rep.empty) break LLoop;
                    client.notify("pushInputs2", rep);
                }
            }
            //Thread.sleep(dur!"msecs"(100));
            while(1){
                auto rep = client.call!(uint[])("pullOutputs2", 128);
                //Thread.sleep(dur!"msecs"(100));
                if(rep.empty) break;
                demod ~= rep;
            }
        }

        writeln(demod);
        cnt += demod.length;
        writefln("Speed: %s", (cast(float)cnt)/((Clock.currTime - start).total!"msecs")*1000*8);
    }
    +/
}
+/
