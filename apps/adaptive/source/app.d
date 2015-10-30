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

import carbon.stream;

import dffdd.filter.diagonal;
import dffdd.filter.lms;
import dffdd.filter.ls;
import dffdd.filter.rls;
import dffdd.filter.mempoly;
import dffdd.filter.polynomial;
import dffdd.utils.fft;

//enum size_t Total = cast(size_t)sampFreq*2;
//enum real sampFreq = 5e6;
enum size_t blockSize = 1024;

//enum size_t N = 4;
//enum size_t P = 3;
//enum size_t Mb = 0;
//enum size_t Mc = 0;
//enum bool withDCBias = true;
//enum bool withIQImbalance = false;

//version = OutputSpectrum;
//version = AdaptiveUseWPH;


void main(string[] args)
{
    JSONValue jv = readText(args[1]).parseJSON();
    dfs(jv, jv["prefix"].str, jv);
}


void dfs(ref JSONValue root, string folder, ref JSONValue jv)
{
    enforce(jv.type == JSON_TYPE.OBJECT);

    foreach(k, ref v; jv.object)
    {
        if(v.type != JSON_TYPE.OBJECT) continue;
        if("ignore" in v.object) continue;

        if("is_dir" in v.object)
            dfs(root, buildPath(folder, k), v);
        else
            implMain(root, folder, k, v);
    }

    return;
}


void implMain(ref JSONValue root, string folder, string fnameId, ref JSONValue jv)
{
    immutable string sendFileNamePrefix = root["send_fn_prefix"].str,
                     recvFileNamePrefix = root["recv_fn_prefix"].str,
                     dataFileNameExt = root["data_fn_ext"].str,
                     sendFileName = buildPath(folder, sendFileNamePrefix ~ fnameId ~ dataFileNameExt),
                     recvFileName = buildPath(folder, recvFileNamePrefix ~ fnameId ~ dataFileNameExt);

    impl(sendFileName, recvFileName, jv["fs"].floating, cast(ptrdiff_t)jv["offset"].integer, cast(size_t)jv["oversampling"].integer);
}


void impl(string sendFileName, string recvFileName, real sampFreq, ptrdiff_t offset, size_t oversampling)
{
    if(!exists(sendFileName) || !exists(recvFileName))
        return;

    File sendFile = File(sendFileName),
         recvFile = File(recvFileName);

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
        //auto adapter = new LMSAdapter!(typeof(state))(state, 0.01, 1024, 0.5);
        auto adapter = makeRLSAdapter(state, 0.999, 1E-7);
        //auto adapter = lsAdapter(state, 8*20);

        return new PolynomialFilter!(typeof(state), typeof(adapter))(state, adapter);
    }();
  }
  else
  {
    auto filter1 = (){
        auto st1 = new PowerState!(cfloat, 8, 1, FilterOptions.usePower)(1),
             st1c = new PowerState!(cfloat, 8, 1, FilterOptions.useConjugate | FilterOptions.usePower)(1),
             st1a = new PowerState!(cfloat, 8, 1, FilterOptions.withConjugate | FilterOptions.usePower)(1),
             st3 = new PowerState!(cfloat, 8, 3, FilterOptions.usePower)(1),
             st3c = new PowerState!(cfloat, 8, 3, FilterOptions.useConjugate | FilterOptions.usePower)(1),
             st3a = new PowerState!(cfloat, 8, 3, FilterOptions.withConjugate | FilterOptions.usePower)(1),
             st5 = new PowerState!(cfloat, 8, 5, FilterOptions.usePower)(1),
             st5c = new PowerState!(cfloat, 8, 5, FilterOptions.useConjugate | FilterOptions.usePower)(1),
             st5a = new PowerState!(cfloat, 8, 5, FilterOptions.withConjugate | FilterOptions.usePower)(1),
             st7 = new PowerState!(cfloat, 8, 7, FilterOptions.usePower)(1),
             st7c = new PowerState!(cfloat, 8, 7, FilterOptions.useConjugate | FilterOptions.usePower)(1),
             st7a = new PowerState!(cfloat, 8, 7, FilterOptions.withConjugate | FilterOptions.usePower)(1),

             st1_2 = new PowerState!(cfloat, 2, 1)(1),
             st1c_2 = new PowerState!(cfloat, 2, 1, FilterOptions.useConjugate)(1),
             st3_2 = new PowerState!(cfloat, 2, 3)(1),
             st3c_2 = new PowerState!(cfloat, 2, 3, FilterOptions.useConjugate)(1),
             //st3_2 = new PowerState!(cfloat, 2, 3, FilterOptions.withConjugate)(1),
             //st5_2 = new PowerState!(cfloat, 2, 5, FilterOptions.withConjugate)(1),
             //st7_2 = new PowerState!(cfloat, 2, 7, FilterOptions.withConjugate)(1),
             //st1c = new PowerState!(cfloat, 8, 1, FilterOptions.useConjugate | FilterOptions.usePower)(1),
             //st3c = new PowerState!(cfloat, 8, 3, FilterOptions.useConjugate | FilterOptions.usePower)(1),
             //st5c = new PowerState!(cfloat, 8, 5, FilterOptions.useConjugate | FilterOptions.usePower)(1),
             //st7c = new PowerState!(cfloat, 8, 7, FilterOptions.useConjugate | FilterOptions.usePower)(1),
             //st7 = new PowerState!(cfloat, 8, 7, FilterOptions.usePower)(1),
             bis = new BiasState!(cfloat)();

        return serialFilter(
                st1_2,
                makeRLSAdapter(st1_2, 0.98, 1E-7),
                //st1c_2,
                //makeRLSAdapter(st1c_2, 0.98, 1E-7),
                //st3_2,
                //makeRLSAdapter(st3_2, 0.999, 1E-7),
                //st3c_2,
                //makeRLSAdapter(st3c_2, 0.999, 1E-7),
                st1,
                //lsAdapter(st1, 400),
                lmsAdapter(st1, 0.0085, 1024, 0.5),
                //makeRLSAdapter(st1, 0.999, 1E-7),
                st1c,
                //lsAdapter(st1c, 1000),
                lmsAdapter(st1c, 0.0085, 1024, 0.5),
                //makeRLSAdapter(st1c, 0.999, 1E-7),
                st3,
                //lsAdapter(st3, 400),
                lmsAdapter(st3, 0.0085, 1024, 0.5),
                //makeRLSAdapter(st3, 0.999, 1E-7),
                st3c,
                //lsAdapter(st3c, 1000),
                lmsAdapter(st3c, 0.0085, 1024, 0.5),
                //makeRLSAdapter(st3c, 0.999, 1E-7),
                st5,
                //lsAdapter(st5, 1000),
                lmsAdapter(st5, 0.0085, 1024, 0.5),
                //makeRLSAdapter(st5, 0.999, 1E-7),
                st5c,
                //lsAdapter(st5c, 1000),
                lmsAdapter(st5c, 0.0085, 1024, 0.5),
                //makeRLSAdapter(st5c, 0.999, 1E-7),
                st7,
                //lsAdapter(st7, 1000),
                lmsAdapter(st7, 0.0085, 1024, 0.5),
                //makeRLSAdapter(st7, 0.999, 1E-7),
                st7c,
                //lsAdapter(st7c, 1000),
                lmsAdapter(st7c, 0.0085, 1024, 0.5),
                //makeRLSAdapter(st7c, 0.999, 1E-7),
                //st7,
                //lsAdapter(st7, 500),
                //lmsAdapter(st7, 1E-1, 1024, 0.5),
                //bis,
                //lmsAdapter(bis, 1E-2, 1024, 0.5),
                );
    }();
  }

    cfloat[] sendBuf = new cfloat[blockSize],
             recvBuf = new cfloat[blockSize],
             intermBuf = new cfloat[blockSize],
             outputBuf = new cfloat[blockSize];

    double[] fftResultRecv = new double[blockSize],
             fftResultSIC = new double[blockSize];

    fftResultRecv[] = 0;
    fftResultSIC[] = 0;

    size_t sumCNT;

    Fft fftObj = new Fft(blockSize);
    auto startTime = Clock.currTime;

    File powerFile = File("errorout_long_%-(%s_%).csv".format(sendFileName.pathSplitter().array()[1 .. $]), "w");

    foreach(blockIdx; -400 .. 400)
    {
        auto sendGets = sendFile.rawRead(sendBuf);
        auto recvGets = recvFile.rawRead(recvBuf);
        assert(sendGets.length == sendBuf.length);
        assert(recvGets.length == recvBuf.length);

        if(blockIdx < 0) continue;

        filter1.apply(sendGets, recvGets, outputBuf);

        if(blockIdx == 300)
        {
            File outFile = File("signalout_%-(%s_%).csv".format(sendFileName.pathSplitter().array()[1 .. $]), "w");
            foreach(i; 0 .. blockSize)
                outFile.writefln("%s,%s,%s,%s,%s,", i, sendGets[i].re, recvGets[i].re, recvGets[i].re - outputBuf[i].re, outputBuf[i].re);
        }

        // fft
        {
            auto outputSpec = fftObj.fftWithSwap(outputBuf);
            auto recvSpec = fftObj.fftWithSwap(recvGets);

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] += (recvSpec[i].re^^2 + recvSpec[i].im^^2);
                fftResultSIC[i] += (outputSpec[i].re^^2 + outputSpec[i].im^^2);
            }
            ++sumCNT;
        }

        if(blockIdx == 0)
        {
            File outFile = File("errorout_start_%-(%s_%).csv".format(sendFileName.pathSplitter().array()[1 .. $]), "w");
            foreach(i; 0 .. blockSize){
                auto pr = recvGets[i].re^^2 + recvGets[i].im^^2,
                     po = outputBuf[i].re^^2 + outputBuf[i].im^^2,
                     c = 10*log10(po / pr);
                outFile.writefln("%s,%s,%s,%s,", i, pr, po, c);
            }
        }
        else if(blockIdx == 300)
        {
            real sumR = 0, sumO = 0;
            File outFile = File("errorout_end_%-(%s_%).csv".format(sendFileName.pathSplitter().array()[1 .. $]), "w");
            foreach(i; 0 .. blockSize){
                auto pr = recvGets[i].re^^2 + recvGets[i].im^^2,
                     po = outputBuf[i].re^^2 + outputBuf[i].im^^2,
                     c = 10*log10(po / pr);

                sumR += pr;
                sumO += po;
                outFile.writefln("%s,%s,%s,%s,", i, pr, po, c);
            }

            writefln("%s, %s [dB]", sendFileName, 10*log10(sumO / sumR));
        }
        
        {
            real sumR = 0, sumO = 0;
            foreach(i; 0 .. blockSize){
                auto pr = recvGets[i].re^^2 + recvGets[i].im^^2,
                     po = outputBuf[i].re^^2 + outputBuf[i].im^^2,
                     c = 10*log10(po / pr);

                sumR += pr;
                sumO += po;
            }

            powerFile.writefln("%s,%s,%s,%s", blockIdx*blockSize, sumR, sumO, 10*log10(sumO / sumR));
        }

        if(blockIdx % 100 == 0){
            real sum = 0; size_t fcnt;
            foreach(i; 0 .. blockSize)
            {
                auto before = 10*log10(fftResultRecv[i] / sumCNT),
                     after = 10*log10(fftResultSIC[i] / sumCNT);

                real freq = i*sampFreq/blockSize-(sampFreq/2);
                //if(abs(freq) < sampFreq / oversampling / 2 * 0.5){
                //    sum += 10.0L ^^(-(before - after)/10);
                //    ++fcnt;
                //}
                if(oversampling == 1 && abs(freq) < sampFreq / 2 * 0.25
                                     && abs(freq) > sampFreq / 10)
                {
                    sum += 10.0L ^^(-(before - after)/10);
                    ++fcnt;
                }else if(abs(freq) < sampFreq / oversampling / 2 * 0.5){
                    sum += 10.0L ^^(-(before - after)/10);
                    ++fcnt;
                }
            }

            if(blockIdx == 300)
            {
                //writefln("%s, %s [dB]", sendFileName, 10*log10(sum / fcnt));

                File outFile = File("fftout_%-(%s_%).csv".format(sendFileName.pathSplitter().array()[1 .. $]), "w");

                foreach(i; 0 .. blockSize)
                {
                    auto before = 10*log10(fftResultRecv[i] / sumCNT),
                         after = 10*log10(fftResultSIC[i] / sumCNT);
                    real freq = i*sampFreq/blockSize-(sampFreq/2);

                    outFile.writefln("%s,%s,%s,%s,%s,%s,", blockIdx, i, freq, before, after, before - after);
                }

                return;
            }

            //writefln("%s [k samples/s],", (blockIdx+1) * blockSize / ((Clock.currTime - startTime).total!"msecs"() / 1000.0L) / 1000.0L);

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] = 0;
                fftResultSIC[i] = 0;
            }
            sumCNT = 0;
        }
    }
}


/+
void main(string[] args)
{
    string sendFilename, recvFilename, outFilename;
    bool bSpeedMode = false, noOutput = false;
    ptrdiff_t seekOffset;

    getopt(args,
        "sendData", &sendFilename,
        "recvData", &recvFilename,
        "outData", &outFilename,
        "speedMode", &bSpeedMode,
        "noOutput", &noOutput,
        "offset", &seekOffset,
    );

    if(outFilename == null && !noOutput){
        auto currTime = Clock.currTime;
        outFilename = format("output_%s_%d_%d_%d_%d_%s_%s_%s.csv", sendFilename.baseName.stripExtension, N, P, Mb, Mc, withDCBias ? "t" : "f", withIQImbalance ? "t" : "f", currTime.toISOString());
    }

    File sendFile = File(sendFilename),
         recvFile = File(recvFilename),
         outFile/* = File(outFilename, "w")*/;

    if(!noOutput)
        outFile = File(outFilename, "w");

    if(seekOffset > 0)
        recvFile.seek(seekOffset * 8);
    else if(seekOffset < 0)
        sendFile.seek(-seekOffset * 8);

    //auto filter1 = {
    //    auto state = new MemoryPolynomialState!(cfloat, N, P, Mb, Mc, withDCBias, withIQImbalance)(1);

    //    //auto adapter = new LSAdapter!(typeof(state))(500);
    //    auto adapter = new LMSAdapter!(typeof(state))(state, 1E-2, 1024, 0.5);

    //    return new PolynomialFilter!(typeof(state), typeof(adapter))(state, adapter);
    //}();
    
    auto st1 = new PowerState!(cfloat, 4, 1, true)(1),
         st3 = new PowerState!(cfloat, 1024, 3, false)(1),
         st5 = new PowerState!(cfloat, 1024, 5, true)(1),
         st7 = new PowerState!(cfloat, 1024, 7, true)(1),
         bis = new BiasState!(cfloat)()
         ;

    auto filter1 = serialFilter(
                                st1,
                                //lsAdapter(st1, 500),
                                lmsAdapter(st1, 1E-2, 1024, 0.5),
                                //st3,
                                //lsAdapter(st3, 500),
                                //lmsAdapter(st3, 1E-3, 1024, 0.5),
                                //st5,
                                //lsAdapter(st5, 500),
                                //lmsAdapter(st5, 1E-3, 1024, 0.5),
                                //st7,
                                //lsAdapter(st7, 500),
                                //lmsAdapter(st7, 1E-3, 1024, 0.5),
                                bis,
                                lmsAdapter(bis, 2E-3, 1024, 0.5),
                                );

    //pragma(msg, filter1.state.state.length * filter1.state.state[0].length);

    cfloat[] sendBuf = new cfloat[blockSize],
             recvBuf = new cfloat[blockSize],
             intermBuf = new cfloat[blockSize],
             outputBuf = new cfloat[blockSize];

    double[] fftResultRecv = new double[blockSize],
             fftResultSIC = new double[blockSize];

    size_t sumCNT;

    Fft fftObj = new Fft(blockSize);
    auto startTime = Clock.currTime;


    foreach(blockIdx; 0 .. Total)
    {
        auto sendGets = sendFile.rawRead(sendBuf);
        auto recvGets = recvFile.rawRead(recvBuf);
        assert(sendGets.length == sendBuf.length);
        assert(recvGets.length == recvBuf.length);

        filter1.apply(sendGets, recvGets, outputBuf);

      version(OutputIteration)
      {
        real sum = 0;
        foreach(i, e; outputBuf){
            sum += e.abs()^^2/blockSize;
            if((blockIdx*blockSize+i) % 10 == 0)
                outFile.writefln("%s,%s,", (blockIdx*blockSize+i)/sampFreq, e.abs()^^2);
        }
        outFile.flush();

        if(blockIdx >= 40) throw new Exception("");
      }

      version(OutputSpectrum)
      {
        // fft
        if(!bSpeedMode){
            auto outputSpec = fftObj.fftWithSwap(outputBuf);
            auto recvSpec = fftObj.fftWithSwap(recvGets);

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] += (recvSpec[i].re^^2 + recvSpec[i].im^^2);
                fftResultSIC[i] += (outputSpec[i].re^^2 + outputSpec[i].im^^2);
            }
            ++sumCNT;
        }

        if(!bSpeedMode && blockIdx % 100 == 0){
            real sum = 0; size_t fcnt;
            foreach(i; 0 .. blockSize)
            {
                auto before = 10*log10(fftResultRecv[i] / sumCNT),
                     after = 10*log10(fftResultSIC[i] / sumCNT);

                real freq = i*sampFreq/blockSize-(sampFreq/2);
                if(abs(freq) < 1500e3){
                    sum += 10.0L ^^(-(before - after)/10);
                    ++fcnt;
                }

                if(!noOutput)
                    outFile.writefln("%s,%s,%s,%s,%s,%s,", blockIdx, i, freq, before, after, before - after);
            }

            writefln("%s,%s,[dB],%s,[k samples/s],", blockIdx, 10*log10(sum / fcnt), (blockIdx+1) * blockSize / ((Clock.currTime - startTime).total!"msecs"() / 1000.0L) / 1000.0L);

            if(!noOutput) outFile.flush();

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] = 0;
                fftResultSIC[i] = 0;
            }
            sumCNT = 0;
        }else if(blockIdx % 1000 == 0)
            writefln("%s,%s,[k samples/s],", blockIdx, (blockIdx+1) * blockSize / ((Clock.currTime - startTime).total!"msecs"() / 1000.0L) / 1000.0L);
      }

        stdout.flush();
    }
}
+/