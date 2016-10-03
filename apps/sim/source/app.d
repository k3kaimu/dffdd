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


auto psdSaveTo(R)(R r, string filename, string resultDir, size_t dropSize, Model model)
{
    alias E = ElementType!R;
    return r.tee(makeSpectrumAnalyzer!E(filename, resultDir, dropSize, model)).toWrappedRange;
}


auto makeSpectrumAnalyzer(C)(string filename, string resultDir, size_t dropSize, Model model)
{
    return makeInstrument!C(delegate void(FiberRange!C r){
        r.drop(dropSize).writePSD(File(buildPath(resultDir, filename), "w"), model.samplingFreq, 1024);
    });
}


auto makeCancellationProbe(C)(string filename, string resultDir, size_t dropSize, Model model)
{
    alias Tup = Tuple!(C, C);

    return makeInstrument!Tup(delegate void(FiberRange!Tup r){
        auto rd = r.drop(dropSize);
        auto ratio = rd.calculateSIC(model.samplingFreq, 1024, model.ofdm.numOfFFT, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling);
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
            auto ratio = rd.calculateSIC(model.samplingFreq, 1024, model.ofdm.numOfFFT, model.ofdm.numOfSubcarrier, model.ofdm.scaleOfUpSampling, 1);
            //if(ratio.isNaN && rd.empty) return;
            if(ratio.isNaN) ratio = 0;

            file.writefln("%s,%s", cnt*1024, 10*log10(ratio));
            file.flush();
            ++cnt;
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

    InputRange!(Complex!real) received;
    {
        auto desired = desiredRandomBits(model).connectToModulator(modOFDM(model), alwaysFalsePointer, model).map!"a*1.0L".psdSaveTo("psd_desired_afterMD.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        if(model.useDTXIQ)  desired = desired.connectToTXIQMixer(model).psdSaveTo("psd_desired_afterTXIQ.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        if(model.useDTXPN)  desired = desired.connectToTXIQPhaseNoise(model).psdSaveTo("psd_desired_afterTXIQPhaseNoise.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        if(model.useDTXPA)  desired = desired.connectToPowerAmplifier(model).psdSaveTo("psd_desired_afterPA.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        if(model.useDTXIQ2) desired = desired.connectToTXIQMixer(model).psdSaveTo("psd_desired_afterTXIQ2.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        desired = desired.connectTo!PowerControlAmplifier((model.thermalNoise.power(model) * (model.SNR + model.lna.NF/* + 4.3 - 10*log10(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling / model.ofdm.numOfSubcarrier)*/).dB.gain^^2).sqrt.V)
                         .psdSaveTo("psd_desired_afterVGA.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);


        auto selfInterference = siRandomBits(model).connectToModulator(modOFDM(model), switchSwapping, model).map!"a*1.0L".psdSaveTo("psd_SI_afterMD.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        if(model.useSTXIQ)  selfInterference = selfInterference.connectToTXIQMixer(model).psdSaveTo("psd_SI_afterTXIQ.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        if(model.useSTXPN)  selfInterference = selfInterference.connectToTXIQPhaseNoise(model).psdSaveTo("psd_SI_afterTXIQPhaseNoise.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        if(model.useSTXPA)  selfInterference = selfInterference.connectToPowerAmplifier(model).psdSaveTo("psd_SI_afterPA.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        if(model.useSTXIQ2) selfInterference = selfInterference.connectToTXIQMixer(model).psdSaveTo("psd_SI_afterTXIQ2.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        selfInterference = selfInterference.drop(model.ofdm.numOfSamplesOf1Symbol/2).connectToMultiPathChannel.connectTo!PowerControlAmplifier((model.thermalNoise.power(model) * (model.INR + model.lna.NF/* + 4.3 - 10*log10(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling / model.ofdm.numOfSubcarrier)*/).dB.gain^^2).sqrt.V)
                                           .psdSaveTo("psd_SI_afterVGA.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);

        received = desired.connectToSwitch(switchDS)
            .add(selfInterference.connectToSwitch(switchSI))
            /*.connectToAWGN(model)*/
            .add(thermalNoise(model).psdSaveTo("psd_thermal_noise.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model))
            .psdSaveTo("psd_rcv_afterAWGN.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);

        if(model.useSRXLN) received = received.connectToLNA(model).psdSaveTo("psd_rcv_afterLNA.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        if(model.useSRXIQ) received = received.connectToRXIQMixer(model).psdSaveTo("psd_rcv_afterRXIQ.csv", resultDir, model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol, model);
        if(model.useSRXQZ) received = received.connectToQuantizer(model);
    }

    auto txReplica = siRandomBits(model).connectToModulator(modOFDM(model), switchSwapping, model).drop(model.ofdm.numOfSamplesOf1Symbol/2);

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
  else static if(filterStructure.endsWith("PH"))
    auto filter = makeParallelHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("CH"))
    auto filter = makeCascadeHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("FHF"))
    auto filter = makeFrequencyHammersteinFilter!filterOptimizer(model);
  else static if(filterStructure.endsWith("WL"))
    auto filter = makeParallelHammersteinFilter!(isOrthogonalized, filterOptimizer, 2)(modOFDM(model), model);
  else static if(filterStructure.endsWith("L"))
    auto filter = makeParallelHammersteinFilter!(isOrthogonalized, filterOptimizer, 1)(modOFDM(model), model);
  else
    static assert("Cannot identify filter model.");
    //auto filter = makeCascadeHammersteinFilter(modOFDM(model), model);
    //auto filter = makeParallelHammersteinFilter(modOFDM(model), model);
    //auto filter = makeParallelHammersteinWithDCMethodFilter(modOFDM(model), model);

    {
        received.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);
        txReplica.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);

      static if(filterStructure.endsWith("FHF"))
      {
        *switchSwapping = true;
        scope(exit) *switchSwapping = false;
      }


        auto fftObj = new Fft(model.blockSize.nextPowOf2(0));

        //File powerFile = File(buildPath(resultDir, "errorout_long.csv"), "w");
        auto sicIterValue = makeCancellationIterationProbe!(Complex!float)("errorout_long.csv", resultDir, 0, model);

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

            // fft
            //{
            //    auto outputSpec = fftObj.fftWithSwap(outps[0 .. $.nextPowOf2(0)]);
            //    auto recvSpec = fftObj.fftWithSwap(recvs[0 .. $.nextPowOf2(0)]);

            //    real sumR = 0, sumO = 0;
            //    foreach(i; 0 .. model.blockSize){
            //        immutable po = outps[i].re^^2 + outps[i].im^^2,
            //                  pr = recvs[i].re^^2 + recvs[i].im^^2;

            //        immutable freq = i * model.samplingFreq / model.blockSize - (model.samplingFreq/2);
            //        if(abs(freq)/model.samplingFreq < (model.ofdm.numOfSubcarrier*1.0)/(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling)
            //        && abs(freq)/model.samplingFreq > 1.0/model.ofdm.numOfFFT/model.ofdm.scaleOfUpSampling){
            //            sumR += pr;
            //            sumO += po;
            //        }
            //    }

            //    powerFile.writefln("%s,%s,%s,%s,", model.blockSize*blockIdx, sumR, sumO, 10*log10(sumO / sumR));
            //}
        }
    }


    {
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
                foreach(i; 0 .. model.blockSize)
                    refrs[i] = recvs[i] - outps[i];

                .put(psdBeforePSD, recvs);
                .put(psdAfterSIC, outps);
                .put(foutPSD, refrs);
                .put(sicValue, recvs.zip(outps));
            }
        }
    }

    received.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);
    txReplica.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);


    if(/*uniform01() > 1*/true)
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

        received.popFrontN(model.numOfPopFront);
        txReplica.popFrontN(model.numOfPopFront);

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

    foreach(blockIdxo; 0 .. 64 * 1024)
    {
        .put(noisePSD, cast(Complex!float)received.front);
        received.popFront();
    }
}


void main()
{
    // ADC&IQ&PA
    foreach(methodName; AliasSeq!("FHF_RLS"/*"FHF", "PH"*//*, "OPH", "OPHDCM", "OCH", "WL", "L",*/ /*"OPHDCM"*/))
        foreach(learningSymbols; [60])
        {
            writeln("START: ", methodName, " : ", learningSymbols);
            writeln("-----------------------------------------");
            scope(exit){
                writeln("-----------------------------------------");
                writeln("END: ", methodName, " : ", learningSymbols);
            }

            Model[] models;
            string[] dirs;

            foreach(inr; /*iota(50, 55, 5)*/[50])
            {
                Model model;
                model.SNR = 5;
                model.INR = inr;
                model.useDTXIQ = false;
                model.useDTXPN = false;
                model.useDTXPA = false;
                model.useSTXIQ = true;
                model.useSTXPN = false;
                model.useSTXPA = true;
                model.useSTXIQ2 = false;
                model.useSRXLN = true;
                model.useSRXIQ = true;
                model.useSRXQZ = true;

              static if(methodName.endsWith("DCM"))
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
                model.swappedSymbols = 100000;
              else
                model.swappedSymbols = 0;

                models ~= model;
                model.numOfPopFront = 0;        // ここ必要？

              static if(methodName.startsWith("FHF"))
                dirs ~= "with_allImpairements_snr%s_inr%s_%s%s_Nswp%s".format(model.SNR, model.INR, methodName, learningSymbols, model.swappedSymbols);
              else
                dirs ~= "with_allImpairements_snr%s_inr%s_%s%s".format(model.SNR, model.INR, methodName, learningSymbols);
            }

            foreach(i; iota(models.length).parallel()){
                mainImpl!methodName(models[i], buildPath("results", dirs[i]));
            }
        }
}
