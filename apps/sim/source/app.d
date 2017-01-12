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
import std.numeric;

import carbon.math : nextPowOf2;

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

import dffdd.mod.primitives;
import dffdd.mod.qam;
import dffdd.mod.ofdm;

import dffdd.dsp.statistics;


import constant;
import models;
import snippet;


real qfunc(real x) { return 0.5 * erfc(x / SQRT2); }


auto psdSaveTo(R)(R r, string filename, string resultDir, size_t dropSize, Model model, Flag!"withRFUpSampling" withRFUpSampling = No.withRFUpSampling)
{
    alias E = ElementType!R;
    return r.tee(makeSpectrumAnalyzer!E(filename, resultDir, dropSize, model, withRFUpSampling)).toWrappedRange;
}


auto makeSpectrumAnalyzer(C)(string filename, string resultDir, size_t dropSize, Model model, Flag!"withRFUpSampling" withRFUpSampling = No.withRFUpSampling)
{
    return makeInstrument!C(delegate void(FiberRange!C r){
        r.drop(dropSize).writePSD(File(buildPath(resultDir, filename), "w"), model.samplingFreq * (withRFUpSampling ? model.scaleOfRFUpSampling : 1), 1024 * (withRFUpSampling ? model.scaleOfRFUpSampling : 1));
    });
}


auto makeCancellationProbe(C)(string filename, string resultDir, size_t dropSize, Model model)
{
    alias Tup = Tuple!(C, C);

    return makeInstrument!Tup(delegate void(FiberRange!Tup r){
        auto rd = r.drop(dropSize);
        auto ratio = rd.calculateSIC(1024, model.ofdm.scaleOfUpSampling);
        auto file = File(buildPath(resultDir, filename), "w");
        file.writefln("%s", 10*log10(ratio));
    });
}


auto makeCancellationIterationProbe(C)(string filename, string resultDir, size_t dropSize, Model model)
{
    alias Tup = Tuple!(C, C);

    return makeInstrument!Tup(delegate void(FiberRange!Tup r){
        auto rd = r.drop(dropSize);
        auto file = File(buildPath(resultDir, filename), "w");

        size_t cnt;
        while(!rd.empty){
            auto ratio = rd.calculateSIC(1024, 1);
            //if(ratio.isNaN && rd.empty) return;
            if(ratio.isNaN) ratio = 0;

            file.writefln("%s,%s", cnt*1024, 10*log10(ratio));
            file.flush();
            ++cnt;
        }
    });
}


auto makeWaveformProbe(C)(string filename, string resultDir, size_t dropSize, Model model)
{
    return makeInstrument!C(delegate void(FiberRange!C r){
        File datfile = File(buildPath(resultDir, filename.stripExtension ~ ".dat"), "w");
        File csvfile = File(buildPath(resultDir, filename.stripExtension ~ ".csv"), "w");
        r = r.drop(dropSize);

        float[2][] chunk = new float[2][](1024);
        while(!r.empty){
            chunk.length = 0;

            foreach(i; 0 .. 1024){
                if(r.empty) break;
                auto frnt = r.front;
                chunk ~= [frnt.re, frnt.im];
                r.popFront();
            }

            datfile.rawWrite(chunk);
            foreach(e; chunk)
                csvfile.writefln("%s,%s,", e[0], e[1]);
        }
    });
}


void mainImpl(string filterType)(Model model, string resultDir)
{
    immutable bool* alwaysFalsePointer = new bool(false);

    immutable ofdmModSignalPower = (){
        auto _modOFDMTest = modOFDM(model);
        return randomBits(1, model).connectToModulator(_modOFDMTest, alwaysFalsePointer, model).measurePower(1024*1024);
    }();

    mkdirRecurse(resultDir);

    bool* switchDS = new bool(false),
          switchSI = new bool(true);

    bool* switchSwapping = new bool(false);

    auto noise = thermalNoise(model);

    immutable numOfModelTrainingSamples = model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol * model.scaleOfRFUpSampling;

    InputRange!(Complex!real) received;
    {
        auto desired = desiredRandomBits(model).connectToModulator(modOFDM(model, Yes.withRFUpSampling), alwaysFalsePointer, model).map!"a*1.0L".psdSaveTo("psd_desired_afterMD.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        if(model.useDTXIQ)  desired = desired.connectToTXIQMixer(model).psdSaveTo("psd_desired_afterTXIQ.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        if(model.useDTXPN)  desired = desired.connectToTXIQPhaseNoise(model).psdSaveTo("psd_desired_afterTXIQPhaseNoise.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        if(model.useDTXPA)  desired = desired.connectToPowerAmplifier(model).psdSaveTo("psd_desired_afterPA.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        if(model.useDTXIQ2) desired = desired.connectToTXIQMixer(model).psdSaveTo("psd_desired_afterTXIQ2.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        desired = desired.connectTo!PowerControlAmplifier((model.thermalNoise.power(model) * (model.SNR + model.lna.NF/* + 4.3 - 10*log10(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling / model.ofdm.numOfSubcarrier)*/).dB.gain^^2).sqrt.V)
                         .psdSaveTo("psd_desired_afterVGA.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);


        auto selfInterference = siRandomBits(model).connectToModulator(modOFDM(model, Yes.withRFUpSampling), switchSwapping, model).map!"a*1.0L".psdSaveTo("psd_SI_afterMD.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        if(model.useSTXIQ)  selfInterference = selfInterference.connectToTXIQMixer(model).psdSaveTo("psd_SI_afterTXIQ.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        if(model.useSTXPN)  selfInterference = selfInterference.connectToTXIQPhaseNoise(model).psdSaveTo("psd_SI_afterTXIQPhaseNoise.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        if(model.useSTXPA)  selfInterference = selfInterference.connectToPowerAmplifier(model).psdSaveTo("psd_SI_afterPA.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        if(model.useSTXIQ2) selfInterference = selfInterference.connectToTXIQMixer(model).psdSaveTo("psd_SI_afterTXIQ2.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        selfInterference = selfInterference/*.drop(model.ofdm.numOfSamplesOf1Symbol * model.scaleOfRFUpSampling/2)*/.connectToMultiPathChannel(model).connectTo!PowerControlAmplifier((model.thermalNoise.power(model) * (model.INR + model.lna.NF/* + 4.3 - 10*log10(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling / model.ofdm.numOfSubcarrier)*/).dB.gain^^2).sqrt.V)
                                           .psdSaveTo("psd_SI_afterVGA.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);

        received = desired.connectToSwitch(switchDS)
            .add(selfInterference.connectToSwitch(switchSI))
            /*.connectToAWGN(model)*/
            .add(thermalNoise(model).psdSaveTo("psd_thermal_noise.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling))
            .psdSaveTo("psd_rcv_afterAWGN.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);

        if(model.useSRXLN) received = received.connectToLNA(model).psdSaveTo("psd_rcv_afterLNA.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        if(model.useSRXIQ) received = received.connectToRXIQMixer(model).psdSaveTo("psd_rcv_afterRXIQ.csv", resultDir, numOfModelTrainingSamples, model, Yes.withRFUpSampling);
        if(model.useSRXQZ) received = received.connectToQuantizer(model);
        received = received.stride(model.scaleOfRFUpSampling).toWrappedRange();
    }

    auto txReplica = siRandomBits(model).connectToModulator(modOFDM(model), switchSwapping, model)/*.drop(model.ofdm.numOfSamplesOf1Symbol/2)*/;

    // モデルの安定化
    *switchDS = false;
    *switchSI = true;
    received.popFrontN(model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol);
    txReplica.popFrontN(model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol);

    // フィルタの学習
    *switchDS = false;
    *switchSI = true;
    auto recvs = new Complex!float[model.blockSize],
         refrs = new Complex!float[model.blockSize],
         outps = new Complex!float[model.blockSize];

    enum string filterStructure = filterType.split("_")[0];
    enum string filterOptimizer = filterType.split("_")[1];

    enum bool isOrthogonalized = filterStructure[0] == 'O';

  static if(filterStructure.endsWith("PHDCM"))
    auto filter = makeParallelHammersteinWithDCMethodFilter!isOrthogonalized(modOFDM(model), model);
  else static if(filterStructure.endsWith("ARPH"))
    auto filter = makeAliasRemovableParallelHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("PH"))
    auto filter = makeParallelHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("CH"))
    auto filter = makeCascadeHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("CWLH"))
    auto filter = makeCascadeWLHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("CWL1H"))
    auto filter = makeCascadeWL1HammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("FHF"))
    auto filter = makeFrequencyHammersteinFilter!filterOptimizer(model);
  else static if(filterStructure.endsWith("WL"))
    auto filter = makeParallelHammersteinFilter!(isOrthogonalized, filterOptimizer, 2)(modOFDM(model), model);
  else static if(filterStructure.endsWith("L"))
    auto filter = makeParallelHammersteinFilter!(isOrthogonalized, filterOptimizer, 1)(modOFDM(model), model);
  else
    static assert("Cannot identify filter model.");


    {
        //received.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);
        //txReplica.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);

      static if(filterStructure.endsWith("FHF"))
      {
        *switchSwapping = true;
        scope(exit) *switchSwapping = false;
      }


        auto fftObj = new Fft(model.blockSize.nextPowOf2(0));

        //File powerFile = File(buildPath(resultDir, "errorout_long.csv"), "w");
        auto sicIterValue = makeCancellationIterationProbe!(Complex!float)("errorout_long.csv", resultDir, 0, model);
        auto waveTransmit = makeWaveformProbe!(Complex!float)("waveform_transmit_init.csv", resultDir, 0, model);
        auto waveBeforeSIC = makeWaveformProbe!(Complex!float)("waveform_beforeSIC_init.csv", resultDir, 0, model);
        auto waveAfterSIC = makeWaveformProbe!(Complex!float)("waveform_afterSIC_init.csv", resultDir, 0, model);
        auto waveFilterOutput = makeWaveformProbe!(Complex!float)("waveform_filterOutput_init.csv", resultDir, 0, model);

        foreach(blockIdx; 0 .. model.numOfFilterTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol / model.blockSize)
        {
          static if(filterStructure.endsWith("FHF"))
          {
            if(blockIdx >= model.swappedSymbols * model.ofdm.numOfSamplesOf1Symbol / model.blockSize)
                *switchSwapping = false;
          }


            foreach(i; 0 .. model.blockSize){
                recvs[i] = (a => complex(a.re, a.im))(received.front);
                refrs[i] = (a => complex(a.re, a.im))(txReplica.front);

                received.popFront();
                txReplica.popFront();
            }

            filter.apply!true(refrs, recvs, outps);
            //filter.apply!false(refrs, recvs, outps);

            if(model.outputWaveform){
                foreach(i; 0 .. model.blockSize){
                    .put(waveTransmit, refrs[i]);
                    .put(waveBeforeSIC, recvs[i]);
                    .put(waveAfterSIC, outps[i]);
                    .put(waveFilterOutput, recvs[i] - outps[i]);
                }
            }

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

            .put(sicIterValue, recvs.zip(outps));
        }
    }


    immutable powerAfterSIC = (){
        real power = 0;

        auto psdBeforePSD = makeSpectrumAnalyzer!(Complex!float)("psd_beforeSIC.csv", resultDir, 0, model);
        auto psdAfterSIC = makeSpectrumAnalyzer!(Complex!float)("psd_afterSIC.csv", resultDir, 0, model);
        auto foutPSD = makeSpectrumAnalyzer!(Complex!float)("psd_filter_output.csv", resultDir, 0, model);
        auto sicValue = makeCancellationProbe!(Complex!float)("cancellation_value.csv", resultDir, 0, model);

        foreach(blockIdxo; 0 .. 1024)
        {
            foreach(i; 0 .. model.blockSize){
                recvs[i] = (a => complex(a.re, a.im))(received.front);
                refrs[i] = (a => complex(a.re, a.im))(txReplica.front);

                received.popFront();
                txReplica.popFront();
            }

            filter.apply!false(refrs, recvs, outps);

            // blockIdxoが10以上になるまで，慣らし運転する
            if(blockIdxo > 10){
                foreach(i; 0 .. model.blockSize){
                    refrs[i] = recvs[i] - outps[i];
                    power += outps[i].re^^2 + outps[i].im^^2;
                }

                .put(psdBeforePSD, recvs);
                .put(psdAfterSIC, outps);
                .put(foutPSD, refrs);
                .put(sicValue, recvs.zip(outps));
            }
        }

        return power / (1024 - 11) / model.blockSize;
    }();

    //received.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);
    //txReplica.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);


    if(model.outputBER)
    {
        // BER count
        *switchDS = true;
        *switchSI = model.withSI;

        size_t inpSize, outSize;
        {
            auto mod = modOFDM(model);
            inpSize = mod.symInputLength;
            outSize = mod.symOutputLength;
        }

        immutable numOfTrainingBits = model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol / outSize * inpSize;
        immutable numOfTrainingSample2 = numOfTrainingBits / inpSize * outSize;

        received.popFrontN(numOfTrainingSample2);
        txReplica.popFrontN(numOfTrainingSample2);

        real berResult = -1;
        auto berCounter = makeInstrument(delegate void (FiberRange!(Complex!float) r){
            auto bits = r
            .connectTo!PowerControlAmplifier(ofdmModSignalPower.sqrt.V)
            .connectToDemodulator(modOFDM(model), model);

            auto refBits = desiredRandomBits(model);
            refBits.popFrontN(numOfTrainingBits);

            // 慣らし運転
            bits.popFrontN(numOfTrainingBits);
            refBits.popFrontN(numOfTrainingBits);

            berResult = measureBER(bits, refBits, model.berCounter.totalBits);
        });

        auto rcvAfterSICPSD = makeSpectrumAnalyzer!(Complex!float)("psd_rcv_afterSIC.csv", resultDir, 0, model);
        auto rcvBeforeSICPSD = makeSpectrumAnalyzer!(Complex!float)("psd_rcv_beforeIC.csv", resultDir, 0, model);

        size_t loopCount;
        while(berResult == -1){
            ++loopCount;

            foreach(i; 0 .. model.blockSize){
                recvs[i] = (a => complex(a.re, a.im))(received.front);
                refrs[i] = (a => complex(a.re, a.im))(txReplica.front);

                received.popFront();
                txReplica.popFront();
            }

            if(model.withSIC && model.withSI)
                filter.apply!false(refrs, recvs, outps);
            else
                outps[] = recvs[];

            foreach(e; outps){
                berCounter.put(Complex!float(e.re, e.im));

                // loopCountが10以上になるまで，慣らし運転する
                if(loopCount > 10)
                    rcvAfterSICPSD.put(e);
            }

            // loopCountが10以上になるまで，慣らし運転する
            if(loopCount > 10)
                foreach(e; recvs)
                    rcvBeforeSICPSD.put(e);
        }

        File file = File(buildPath(resultDir, "ber.csv"), "w");
        file.writeln(berResult);

        writefln("%s,%s,%2.0f,%2.0f,%s", model.withSI ? 1 : 0, model.withSIC ? 1 : 0, model.SNR, model.INR, berResult);
    }

    // ノイズ電力
    auto noisePSD = makeSpectrumAnalyzer!(Complex!float)("psd_noise_floor.csv", resultDir, 0, model);

    //model.lna.GAIN = 0;
    //received = thermalNoise(model).connectToLNA(model).toWrappedRange;

    *switchSI = false;
    *switchDS = false;
    received.popFrontN(1024);

    real noisePower = 0;
    foreach(blockIdxo; 0 .. 64 * 1024)
    {
        auto c = received.front;
        .put(noisePSD, cast(Complex!float)c);
        noisePower += c.re^^2 + c.im^^2;
        received.popFront();
    }
    noisePower /= 64 * 1024;

    File(buildPath(resultDir, "remained-interference-to-noise-power-ratio.csv"), "w")
    .writeln(10*log10((powerAfterSIC - noisePower)/(noisePower)));
}


void mainJob()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    auto taskList = new MultiTaskList();

    // ADC&IQ&PA
    foreach(methodName; AliasSeq!("FHF_LS", "PH_LS"))
        foreach(learningSymbols; [60])
        {
            Model[] models;
            string[] dirs;

            foreach(inr; iota(20, 85, 5))
            foreach(iip3; [17])
            // foreach(irr; [15, 25, 35])
            // foreach(txp; iota(10, 32, 2))
            // foreach(bUseIQ2; [false, true])
            {
                immutable bUseIQ2 = true;
                Model model;
                model.SNR = 20;
                model.INR = inr;
                // model.pa.TX_POWER = txp;
                // model.txIQMixer.IIR = irr;
                // model.rxIQMixer.IIR = irr;
                model.quantizer.numOfBits = 14;
                model.pa.IIP3 = iip3;

                // 再現する非線形性の選択
                model.useDTXIQ = false;
                model.useDTXPN = false;
                model.useDTXPA = false;
                model.useSTXIQ = !bUseIQ2;
                model.useSTXPN = false;
                model.useSTXPA = true;
                model.useSTXIQ2 = bUseIQ2;
                model.useSRXLN = true;
                model.useSRXIQ = true;
                model.useSRXQZ = true;

                // ベースバンド信号波形の出力
                model.outputWaveform = false;

                model.channel.taps = 64;
                model.firFilter.taps = 64;

                model.outputBER = false;

              static if(methodName.split("_")[0].endsWith("DCM"))
              {
                if(learningSymbols == 10){
                    model.learningSymbols = 4;
                    model.learningCount = 3;
                }else{
                    model.learningSymbols = 4;
                    model.learningCount = 15;
                }
              }
              else
              {
                model.learningSymbols = learningSymbols;
                model.learningCount = 1;
              }


              static if(methodName.startsWith("FHF"))
              {
                model.swappedSymbols = 100000;
                model.numOfFilterTrainingSymbols = 100;
              }
              else
              {
                model.swappedSymbols = 0;
                model.numOfFilterTrainingSymbols = 1000;
              }

              static if(methodName.endsWith("LMS") || methodName.endsWith("RLS"))
                model.numOfFilterTrainingSymbols = 100;

                model.ofdm.scaleOfUpSampling = 4;

                models ~= model;

              static if(methodName.startsWith("FHF"))
                dirs ~= "IIP3_%s_inr%s_%s%s_Nswp%s".format(iip3, model.INR, methodName, learningSymbols, model.swappedSymbols);
              else
                dirs ~= "IIP3_%s_inr%s_os%s_%s%s".format(iip3, model.INR, model.ofdm.scaleOfUpSampling, methodName, learningSymbols);
            }

            foreach(i; 0 .. models.length)
                taskList.append((Model m, string dir){ mainImpl!methodName(m, dir); }, models[i], buildPath("results", dirs[i]));
        }

    //auto scheduler = new MPITaskScheduler();
    jobRun(taskList);
}


void main()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    //jobRun(1, 0, {
        mainJob();
    //});
}
