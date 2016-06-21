import std.algorithm;
import std.array;
import std.complex;
import std.datetime;
import std.format;
import std.getopt;
import std.math;
import std.numeric;
import std.path;
import std.range;
import std.stdio;
import std.stdio;
import std.string;
import std.json;
import std.file;
import std.exception;
import std.typetuple;

import carbon.stream;

import dffdd.filter.diagonal;
import dffdd.filter.lms;
import dffdd.filter.ls;
import dffdd.filter.rls;
import dffdd.filter.mempoly;
import dffdd.filter.polynomial;
import dffdd.utils.fft;
import dffdd.filter.orthogonalize;
import dffdd.utils.jsvisitor;

//enum size_t Total = cast(size_t)sampFreq*2;
//enum real sampFreq = 5e6;
enum ptrdiff_t blockSize = 1024;

//enum size_t N = 4;
//enum size_t P = 3;
//enum size_t Mb = 0;
//enum size_t Mc = 0;
//enum bool withDCBias = true;
//enum bool withIQImbalance = false;

//version = OutputSpectrum;
//version = AdaptiveUseWPH;
//version = OutputCancelledSignal;

alias BasisFunctions = TypeTuple!(x => x,
                                  x => x.conj,
                                  x => x * (x.re^^2 + x.im^^2),
                                  x => x.conj * (x.re^^2 + x.im^^2),
                                  x => x * (x.re^^2 + x.im^^2)^^2,
                                  x => x.conj * (x.re^^2 + x.im^^2)^^2,
                                  x => x * (x.re^^2 + x.im^^2)^^4,
                                  x => x.conj * (x.re^^2 + x.im^^2)^^4,
                                  );


void main(string[] args)
{
    static struct Visitor
    {
        string resultDir;


        Visitor onDirectory(VisitorData data)
        {
            return this;
        }


        void onFile(VisitorData data)
        {
            auto dir = this.resultDir.buildPath(buildPath(data.parentKeys));
            impl(dir, data.txFileName, data.rxFileName, data.samplingFreq, data.offsetTX, data.oversampling, data.parentKeys.length);
        }
    }

    Visitor v;
    v.resultDir = expandTilde("~/exp_result").buildPath((x => [x[2][2..$],x[0],x[1],x[3],x[4],x[5]].join("_"))((x => x.filter!"a.length".array)(__DATE__.split(" ") ~ __TIME__.split(":"))));

    auto rootJSON = args[1].readText().parseJSON();
    visitJSON(rootJSON, v);
}


void impl(string resultDir, string sendFileName, string recvFileName, real sampFreq, ptrdiff_t offset, size_t oversampling, size_t depth)
{
    //writeln(sendFileName);
    //writeln(recvFileName);
    if(!exists(sendFileName) || !exists(recvFileName))
        return;

    mkdirRecurse(resultDir);

    File sendFile = File(sendFileName),
         recvFile = File(recvFileName);

    // offsetのところで相関ピークだけど，一個後ろのサンプルにも相関値がいくらかあるはず．
    offset -= 1;    // txからみて-1, rxからみて+1
    //offset += 1;

    if(offset > 0)
        recvFile.seek(offset * 8);
    else if(offset < 0)
        sendFile.seek(-offset * 8);

    //sendFile.seek(1_000_000 * 8);
    //recvFile.seek(1_000_000 * 8);


  version(AdaptiveUseWPH)
  {
    auto filter1 = {
        auto state = new MemoryPolynomialState!(cfloat, 8, 4, 0, 0, false, true)(1);

        //writeln("Use");
        //auto adapter = new LMSAdapter!(typeof(state))(state, 0.001, 1024, 0.5);
        auto adapter = makeRLSAdapter(state, 1, 1E-7);
        //auto adapter = lsAdapter(state, 8*20);

        return new PolynomialFilter!(typeof(state), typeof(adapter))(state, adapter);
    }();
  }
  else
  {

  //version(AdaptiveUseWPH)
    //enum filterName = "parallelFilter";
  //else
    //enum filterName = "serialFilter";


    auto orthogonalizer = new GramSchmidtOBFFactory!(cfloat, BasisFunctions)();

    {
        orthogonalizer.start();
        File sendFile2 = File(sendFileName);
        sendFile2.seek(1024*8*100);
        scope(exit)
            orthogonalizer.finish();

        cfloat[] buf = new cfloat[1024];

        foreach(i; -400 .. 400){
            ptrdiff_t n = sendFile2.rawRead(buf).length;
            if(i < 0) continue;

            orthogonalizer.put(buf[0 .. n]);
        }
    }

    cfloat[][BasisFunctions.length] coefs;
    foreach(i, ref e; coefs){
        e = new cfloat[BasisFunctions.length];
        orthogonalizer.getCoefs(i, e);
        //foreach(j, ref ee; e){
        //    if(i == j) ee = 1+0i;
        //    else ee = 0+0i;
        //}
    }

    //writeln(coefs);
    //if(coefs[0]) return;
    //coefs[0][0] = 1+0i;
    //writeln(coefs[0]);

    auto filter1 = (){
        auto st1 = (new FIRFilter!(cfloat, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[0].dup);
        auto st12 = (new FIRFilter!(cfloat, 3, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[0].dup);
        auto st1c = (new FIRFilter!(cfloat, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[1].dup);
        //auto st2 = (new FIRFilter!(cfloat, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[2].dup);
        auto st3 = (new FIRFilter!(cfloat, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[2].dup);
        //auto st32 = (new FIRFilter!(cfloat, 3, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[2].dup);
        auto st3c = (new FIRFilter!(cfloat, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[3].dup);
        auto st5 = (new FIRFilter!(cfloat, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[4].dup);
        auto st5c = (new FIRFilter!(cfloat, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[5].dup);
        auto st7 = (new FIRFilter!(cfloat, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[6].dup);
        auto st7c = (new FIRFilter!(cfloat, 8, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[7].dup);
        //auto st5 = (new FIRFilter!(cfloat, 16, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[4].dup);
        //auto st7 = (new FIRFilter!(cfloat, 16, true)).inputTransformer!((x, h) => OBFEval!BasisFunctions(x, h))(coefs[5].dup);

        return serialFilter(
                //st12,
                //makeRLSAdapter(st12, 0.98, 1E-7),
                //st32,
                //makeRLSAdapter(st32, 0.98, 1E-7),
                //st1c_2,
                //makeRLSAdapter(st1c_2, 0.98, 1E-7),
                //st3_2,
                //makeRLSAdapter(st3_2, 0.92, 1E-7),
                //st3c_2,
                //makeRLSAdapter(st3c_2, 0.999, 1E-7),
                st1,
                //lsAdapter(st1, 50000),
                lmsAdapter(st1, 0.001, 1024, 0.5),
                //makeRLSAdapter(st1, 0.99, 1E-7),
                st1c,
                //lsAdapter(st1c, 1000),
                lmsAdapter(st1c, 0.001, 1024, 0.5),
                //makeRLSAdapter(st1c, 1 - 1E-4, 1E-7),
                //st2,
                //lsAdapter(st2, 5000),
                //lmsAdapter(st2, 0.002, 1024, 0.5),
                //makeRLSAdapter(st2, /*0.9997*/1, 1E-7),
                st3,
                //lsAdapter(st3, 50000),
                lmsAdapter(st3, 0.001, 1024, 0.5),
                //makeRLSAdapter(st3, 1 - 1E-4, 1E-7),
                st3c,
                //lsAdapter(st3c, 1000),
                lmsAdapter(st3c, 0.001, 1024, 0.5),
                //makeRLSAdapter(st3c, 1 - 1E-4, 1E-7),
                //st4,
                //lsAdapter(st4, 5000),
                //lmsAdapter(st4, 0.002, 1024, 0.5),
                //makeRLSAdapter(st4, /*0.9997*/1, 1E-7),
                //st4a,
                //lmsAdapter(st4a, 0.0005, 1024, 0.5),
                st5,
                //lsAdapter(st5, 5000),
                lmsAdapter(st5, 0.001, 1024, 0.5),
                //makeRLSAdapter(st5, 1 - 1E-4, 1E-7),
                st5c,
                //lsAdapter(st5c, 1000),
                lmsAdapter(st5c, 0.001, 1024, 0.5),
                //makeRLSAdapter(st5c, 1 - 1E-4, 1E-7),
                st7,
                //lsAdapter(st7, 5000),
                lmsAdapter(st7, 0.001, 1024, 0.5),
                //makeRLSAdapter(st7, 1 - 1E-6, 1E-7),
                //st12,
                //lsAdapter(st1, 5000),
                //lmsAdapter(st12, 0.010, 1024, 0.5),
                //makeRLSAdapter(st7, 0.9997, 1E-7),
                st7c,
                //lsAdapter(st7c, 1000),
                lmsAdapter(st7c, 0.001, 1024, 0.5),
                //makeRLSAdapter(st7c, 1 - 1E-6, 1E-7),
                //st5,
                //lmsAdapter(st5, 0.0085, 1024, 0.5),
                //makeRLSAdapter(st5, 1, 1E-7),
                //st3,
                //lmsAdapter(st3, 0.01, 1024, 0.5),
                //makeRLSAdapter(st3, 1, 1E-7),
                //st1,
                //lmsAdapter(st1, 0.01, 1024, 0.5),
                //makeRLSAdapter(st1, 1, 1E-7),
                );
    }();
  }

    cfloat[] sendBuf = new cfloat[blockSize],
             recvBuf = new cfloat[blockSize],
             intermBuf = new cfloat[blockSize],
             outputBuf = new cfloat[blockSize],
             interBuf1 = new cfloat[blockSize],
             interBuf2 = new cfloat[blockSize];

    double[] fftResultRecv = new double[blockSize],
             fftResultSIC = new double[blockSize],
             fftResultFilterOutput = new double[blockSize];

    fftResultRecv[] = 0;
    fftResultSIC[] = 0;

    size_t sumCNT;

    Fft fftObj = new Fft(blockSize);
    auto startTime = Clock.currTime;

    //immutable string fileSuffix = format("%-(%s_%)", sendFileName.stripExtension.pathSplitter.array[$ >= depth ? $-depth : 0 .. $]);

    File powerFile = File(buildPath(resultDir, "errorout_long.csv"), "w");
    version(OutputCancelledSignal) File cancelledFile = File(buildPath(resultDir, "cancelled_signal.dat"), "w");

  version(OutputCancelledSignal)
  {
    immutable ptrdiff_t startIdx = -400*1024/blockSize,
                        endIdx = ptrdiff_t.max,
                        learningEndIdx = 400*1024/blockSize;
  }
  else
  {
    immutable ptrdiff_t startIdx = -400*1024/blockSize,
                        endIdx = 400*1024/blockSize,
                        learningEndIdx = 400*1024/blockSize;
  }

    foreach(blockIdx; startIdx .. endIdx)
    {
        auto sendGets = sendFile.rawRead(sendBuf);
        auto recvGets = recvFile.rawRead(recvBuf);
        //enforce(sendGets.length == sendBuf.length);
        if(sendGets.length != sendBuf.length) return; // end
        //enforce(recvGets.length == recvBuf.length);
        if(recvGets.length != recvBuf.length) return; // end

        if(blockIdx < 0) continue;

      version(OutputCancelledSignal)
      {
        if(blockIdx % (100*1024/blockSize) == 0)
            writefln("idx=%s, %s [k samples/s],", blockIdx, (blockIdx+1) * blockSize / ((Clock.currTime - startTime).total!"msecs"() / 1000.0L) / 1000.0L);
      }

        //writeln(blockIdx);
        //filter7.apply(sendGets, recvGets, interBuf1);
        //filter5.apply(sendGets, interBuf1, interBuf2);
        //filter3.apply(sendGets, interBuf2, interBuf1);
        if(blockIdx < learningEndIdx)
          filter1.apply!true(sendGets, recvGets, outputBuf);
        else
          filter1.apply!false(sendGets, recvGets, outputBuf);

        if(blockIdx == 300*1024/blockSize)
        {
            File outFile = File(buildPath(resultDir, "signalout.csv"), "w");
            foreach(i; 0 .. blockSize)
                outFile.writefln("%s,%s,%s,%s,%s,", i, sendGets[i].re, recvGets[i].re, recvGets[i].re - outputBuf[i].re, outputBuf[i].re);
        }

        // fft
        {
            auto outputSpec = fftObj.fftWithSwap(outputBuf);
            auto recvSpec = fftObj.fftWithSwap(recvGets);

            real sumR = 0, sumO = 0;
            foreach(i; 0 .. blockSize){
                immutable real pr = recvSpec[i].re^^2 + recvSpec[i].im^^2,
                               po = outputSpec[i].re^^2 + outputSpec[i].im^^2;

                fftResultRecv[i] += pr;
                fftResultSIC[i] += po;
                fftResultFilterOutput[i] += (x => x.re^^2 + x.im^^2)(recvSpec[i] - outputSpec[i]);
                real freq = i*sampFreq/blockSize-(sampFreq/2);

                if(abs(freq) <= (sampFreq / 2 / oversampling * 0.80)){
                    if(abs(freq) >= (sampFreq / 64)){
                        sumR += pr;
                        sumO += po;
                    }
                }
            }

            powerFile.writefln("%s,%s,%s,%s,", blockIdx*blockSize, sumR, sumO, 10*log10(sumO / sumR));
            ++sumCNT;
        }

        /*
        波形出力の場合
        */
        version(OutputCancelledSignal){
            // 400以上は波形出力したら終わる
            if(blockIdx > 400*1024/blockSize){
                cancelledFile.rawWrite(outputBuf);
                continue;
            }
        }

        if(blockIdx == 0)
        {
            File outFile = File(buildPath(resultDir, "errorout_start.csv"), "w");
            foreach(i; 0 .. blockSize){
                if(i > 20000) break;
                auto pr = recvGets[i].re^^2 + recvGets[i].im^^2,
                     po = outputBuf[i].re^^2 + outputBuf[i].im^^2,
                     c = -1*10*log10(po / pr);

                outFile.writefln("%s,%s,%s,%s,", i, pr, po, c);
            }
        }

        if(blockIdx % (100*1024/blockSize) == 0){
            if(blockIdx == 300*1024/blockSize)
            {
                File outFile = File(buildPath(resultDir, "fftout.csv"), "w");

                real sumR = 0, sumO = 0;
                foreach(i; 0 .. blockSize)
                {
                    auto pr = fftResultRecv[i],
                         po = fftResultSIC[i],
                         before = 10*log10(fftResultRecv[i] / sumCNT),
                         after = 10*log10(fftResultSIC[i] / sumCNT),
                         filterOutput = 10*log10(fftResultFilterOutput[i] / sumCNT);
                    real freq = i*sampFreq/blockSize-(sampFreq/2);

                    if(abs(freq) <= (sampFreq / 2 / oversampling * 0.80)){
                        if(abs(freq) >= (sampFreq / 64)){
                            sumR += pr;
                            sumO += po;
                        }
                    }

                    outFile.writefln("%s,%s,%s,%s,%s,%s,%s", blockIdx, i, freq, filterOutput, before, after, before - after);
                }

                writefln("%s, %s [dB](%s / %s = %s)", sendFileName, 10*log10(sumO / sumR), sumR, sumO, sumR / sumO);
            }

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] = 0;
                fftResultSIC[i] = 0;
                fftResultFilterOutput[i] = 0;
            }
            sumCNT = 0;
        }
    }
}
