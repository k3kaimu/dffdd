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
        selfInterference = selfInterference/*.drop(model.ofdm.numOfSamplesOf1Symbol/2)*/.connectToMultiPathChannel(model).connectTo!PowerControlAmplifier((model.thermalNoise.power(model) * (model.INR + model.lna.NF/* + 4.3 - 10*log10(model.ofdm.numOfFFT * model.ofdm.scaleOfUpSampling / model.ofdm.numOfSubcarrier)*/).dB.gain^^2).sqrt.V)
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

    auto txReplica = siRandomBits(model).connectToModulator(modOFDM(model), switchSwapping, model)/*.drop(model.ofdm.numOfSamplesOf1Symbol/2)*/;
    auto orthTrainReplica = siRandomBits(model).connectToModulator(modOFDM(model), switchSwapping, model);

    // モデルの安定化
    *switchDS = false;
    *switchSI = true;
    received.popFrontN(model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol);
    txReplica.popFrontN(model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol);
    orthTrainReplica.popFrontN(model.numOfModelTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol);

    // フィルタの学習
    *switchDS = false;
    *switchSI = true;
    auto recvs = new Complex!float[model.blockSize],
         refrs = new Complex!float[model.blockSize],
         outps = new Complex!float[model.blockSize];

  enum string filterStructure = filterType.split("_")[0];
  enum string filterOptimizer = filterType.split("_")[1];

//   enum bool isOrthogonalized = filterStructure[0] == 'O';

  static if(filterStructure.endsWith("PHDCM"))
    auto filter = makeParallelHammersteinWithDCMethodFilter!isOrthogonalized(modOFDM(model), model);
  else static if(filterStructure.endsWith("ARPH"))
    auto filter = makeAliasRemovableParallelHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("PH"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("CH"))
    auto filter = makeCascadeHammersteinFilter!(filterOptimizer)(modOFDM(model), model);
//   else static if(filterStructure.endsWith("CWLH"))
//     auto filter = makeCascadeWLHammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
//   else static if(filterStructure.endsWith("CWL1H"))
//     auto filter = makeCascadeWL1HammersteinFilter!(isOrthogonalized, filterOptimizer)(modOFDM(model), model);
  else static if(filterStructure.endsWith("CFHF"))
    auto filter = makeFrequencyCascadeHammersteinFilter!(filterOptimizer)(model);
  else static if(filterStructure.endsWith("FHF"))
    auto filter = makeFrequencyHammersteinFilter!filterOptimizer(model);
  else static if(filterStructure.endsWith("WL"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, 2)(modOFDM(model), model);
  else static if(filterStructure.endsWith("L"))
    auto filter = makeParallelHammersteinFilter!(filterOptimizer, 1)(modOFDM(model), model);
  else
    static assert("Cannot identify filter model.");

//   static if(is(typeof(filter.learningFromTX([Complex!float.init]))))
  {
    if(model.orthogonalizer.enabled && model.orthogonalizer.numOfTrainingSymbols != 0)
    {
        static if(filterStructure.endsWith("FHF")){
            writeln("SWAP");
            *switchSwapping = true;
        }

        // 学習系列の生成
        Complex!float[] txs;
        foreach(i; 0 .. model.orthogonalizer.numOfTrainingSymbols * model.ofdm.numOfSamplesOf1Symbol){
            static if(filterStructure.endsWith("FHF"))
            if(i >= model.swappedSymbols * model.ofdm.numOfSamplesOf1Symbol){
                if(*switchSwapping) writeln("END SWAP");
                *switchSwapping = false;
            }

            txs ~= orthTrainReplica.front;
            orthTrainReplica.popFront();
        }

        // 学習させる
        filter.learningFromTX(txs);
        writeln("Training of the distorter is finished.");
    }
  }

    {
        //received.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);
        //txReplica.popFrontN(model.ofdm.numOfSamplesOf1Symbol/2*5);

      static if(filterStructure.endsWith("FHF"))
      {
        *switchSwapping = true;
        //scope(exit) *switchSwapping = false;
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

            filter.apply!(Yes.learning)(refrs, recvs, outps);
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

            filter.apply!(No.learning)(refrs, recvs, outps);

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
                filter.apply!(No.learning)(refrs, recvs, outps);
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


void mainJob()
{
    //import tuthpc.mpi;
    import tuthpc.taskqueue;

    auto taskList = new MultiTaskList();

    // ADC&IQ&PA
    foreach(methodName; AliasSeq!(/*"FHF_LS",*/ /*"OFHF_LS",*/ "CFHF_LS",
                // "FHF_LMS", "FHF_LS", "OPH_LS", "OPH_RLS", "OPH_LMS", "OCH_LS", "OCH_RLS", "OCH_LMS", "WL_LS", "WL_RLS", "WL_LMS", "L_LS", "L_RLS", "L_LMS" /*"FHF", "PH"*//*, "OPH", "OPHDCM", "OCH", "WL", "L",*/ /*"OPHDCM"*/
            ))
        foreach(learningSymbols; iota(60, 65, 5)) foreach(orthTrainingSymbols; [10000])
        {
            Model[] models;
            string[] dirs;

            foreach(inr; iota(50, 55, 5)) foreach(txp; iota(15, 18, 3))
            {
                Model model;
                model.SNR = 20;
                model.INR = inr;
                model.pa.TX_POWER = txp;

                // 再現する非線形性の選択
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

                model.orthogonalizer.numOfTrainingSymbols = orthTrainingSymbols;

                if(methodName[0] == 'O')
                    model.orthogonalizer.enabled = true;
                else
                    model.orthogonalizer.enabled = false;

                // ベースバンド信号波形の出力
                model.outputWaveform = false;

                model.channel.taps = 1;
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


              static if(methodName.split("_")[0].endsWith("FHF"))
              {
                model.swappedSymbols = 100000;
                model.numOfFilterTrainingSymbols = 100;
              }
              else
              {
                model.swappedSymbols = 0;
                model.numOfFilterTrainingSymbols = 1000;
              }

                models ~= model;

              static if(methodName.endsWith("LMS") || methodName.endsWith("RLS"))
                model.numOfFilterTrainingSymbols = 100;

              static if(methodName.split("_")[0].endsWith("FHF"))
                dirs ~= "TXP%s_inr%s_%s%s_Nswp%s_orth%s".format(model.pa.TX_POWER, model.INR, methodName, learningSymbols, model.swappedSymbols, model.orthogonalizer.numOfTrainingSymbols);
              else
                dirs ~= "TXP%s_inr%s_%s%s_orth%s".format(model.pa.TX_POWER, model.INR, methodName, learningSymbols, model.orthogonalizer.numOfTrainingSymbols);
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
