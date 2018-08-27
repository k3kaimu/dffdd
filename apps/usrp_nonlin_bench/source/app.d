module app;

import std.stdio;

import carbon.channel;
import core.thread;

import std.algorithm;
import std.bitmanip;
import std.complex;
import std.concurrency;
import std.conv;
import std.datetime;
import std.exception;
import std.math;
import std.random;
import std.range;
import std.string;
import std.meta;

import uhd.capi.types.metadata;
import uhd.usrp;

import dffdd.blockdiagram.noise;
import dffdd.dsp.convolution;
import dffdd.mod.ofdm;
import dffdd.mod.qpsk;
import dffdd.mod.bpsk;
import dffdd.utils.fft;

//import zmqd;
import preamble;
import spec;
import impl;

import core.memory;

shared bool endFlag;

enum Mode { send, recv }
enum ModMethod { bpsk, qpsk }


enum size_t unitBits = 1024;
enum size_t totalBits = unitBits*10_000;
enum bool loadFromFile = false;


void main(string[] args)
{
    immutable usrpAddr = args[1];
    immutable clockSource = args[2];
    immutable sendOrReceive = args[3] == "send" ? Mode.send : Mode.recv;
    immutable modmethod = args[4] == "bpsk" ? ModMethod.bpsk : ModMethod.qpsk;
    immutable paygain_dB = to!real(args[5]);
    immutable paygain_dB_str = args[5];
    immutable isSimulation = (args.length > 6) && args[6] == "sim";
    immutable batchSize = isSimulation ? 1000 : 100;   // # of transmit bits = 1024 * batchSize

    immutable paygain = 10.0L^^(paygain_dB / 20);

    Random rnd;
    rnd.seed(0xDEADBEEF);
    immutable data = () {
        return repeat(&rnd).take((unitBits * batchSize) / (size_t.sizeof * 8)).map!(rnd => uniform!size_t(*rnd)).array().assumeUnique;
    }();

  static if(!loadFromFile)
    foreach(i, Mod; AliasSeq!(BPSK!(Complex!float), QPSK!(Complex!float))){
        if(modmethod.to!string == Mod.stringof[0..4].toLower){
            Mod mod;

            if(isSimulation)
            {
                auto channel = channel!ICA();

                auto zmqTX = spawn(
                                &txMain!Mod,
                                channel,
                                data,
                                mod,
                                paygain,
                                isSimulation);

                auto zmqRX = spawn(
                                &rxMain!Mod,
                                channel,
                                data,
                                mod,
                                paygain_dB_str,
                                isSimulation);
            }
            else
            {
                auto txChannel = channel!ICA();
                auto rxChannel = channel!ICA();

                if(sendOrReceive == Mode.send){
                    auto txtid = spawn(
                                    &transmitWorker,
                                    usrpAddr,
                                    (100e6/2/25),
                                    15,
                                    5.015e9,
                                    clockSource,
                                    txChannel);

                    auto zmqTX = spawn(
                                    &txMain!Mod,
                                    txChannel,
                                    data,
                                    mod,
                                    paygain,
                                    isSimulation);
                }else{
                    auto rxtid = spawn(
                                    &receiveWorker,
                                    usrpAddr,
                                    (100e6/2/25),
                                    40,
                                    5.015e9,
                                    clockSource,
                                    rxChannel);

                    auto zmqRX = spawn(
                                    &rxMain!Mod,
                                    rxChannel,
                                    data,
                                    mod,
                                    paygain_dB_str,
                                    isSimulation);
                }
            }
        }
    }


    static if(loadFromFile)
    {
        alias Mod = BPSK;
        BPSK mod;

        immutable signal = (){
                float[] filedata = new float[1024];
                File file = File("AWGN_yamada/%.1fdB_AWGN.bin".format(round(10*log10(snr))), "r");
                file.rawRead(filedata);

                File csv = File("data.csv", "w");
                foreach(e; filedata) csv.writefln("%s", e);

                return filedata.map!(a => Complex!float(a, 0)).repeat(batchSize).joiner.array().assumeUnique;
            }();

        if(signal !is null) return;

        writeln("LOAD FROM FILE: ", signal.length);

        if(sendOrReceive == Mode.send){
            auto txChannel = channel!ICA();

            auto txtid = spawn(
                            &transmitWorker,
                            usrpAddr,
                            (100e6/2/25),
                            15,
                            5.015e9,
                            clockSource,
                            txChannel);

            auto zmqTX = spawn(
                            &txMain_loadFromFile!Mod,
                            txChannel,
                            signal);
        }else{
            auto rxChannel = channel!ICA();

            auto rxtid = spawn(
                            &receiveWorker,
                            usrpAddr,
                            (100e6/2/25),
                            40,
                            5.015e9,
                            clockSource,
                            rxChannel);

            auto zmqRX = spawn(
                            &rxMain_loadFromFile!Mod,
                            rxChannel,
                            signal,
                            snr);
        }
    }


    foreach(line; stdin.byLine){
        if(line.startsWith("s")){
        }else if(line.startsWith("e")){
            endFlag = true;
            break;
        }
    }
}


void transmitWorker(
    string addr,
    double rate,
    double gain,
    double freq,
    string clockSource,
    shared(Channel!ICA) channel)
{
    try{
        auto usrp = USRP(MultiDeviceAddress([addr]));

        usrp.txRate = rate;
        writeln("Actual TX Rate: ", usrp.txRate);

        usrp.txGain = gain;
        writeln("Actual TX Gain: ", usrp.txGain);

        usrp.txFreq = freq;
        writefln("Actual TX freq: %s [MHz]", usrp.txFreq / 1e6);

        usrp.clockSource = clockSource;
        writefln("Actual clock source: %s", usrp.clockSource);

        // usrp.timeSource = clockSource;
        // writefln("Actual clock source: %s", usrp.timeSource);

        auto txStreamer = usrp.makeTxStreamer(StreamArgs("fc32", "sc16", "", [0]));
        auto txMaxlen = txStreamer.maxNumSamps;
        writefln("TXBuffer size in samples: %s", txMaxlen);

        auto txmdstart = TxMetaData(false, 0, 0.1, true, false);
        auto txmdmiddle = TxMetaData(false, 0, 0.1, false, false);
        auto txmdend = TxMetaData(false, 0, 0.1, false, true);

        bool isSendMode = false;

        while(!endFlag){
            if(auto p = channel.pop!ICA){
                writeln("POP: ", (*p).length);

                if(!isSendMode && (*p).length == 0){
                    isSendMode = true;
                    //txStreamer.send(null, txmdstart, 0.00001);
                    continue;
                }else if(isSendMode && (*p).length == 0){
                    isSendMode = false;
                    //txStreamer.send(null, txmdend, 0.00001);
                    continue;
                }else if(isSendMode && (*p).length != 0){
                    auto chunk = *p;
                    writeln("START");
                    txStreamer.send(null, txmdstart, 0.0001);

                    size_t cnt;
                    while(cnt < chunk.length)
                        cnt += txStreamer.send(chunk[cnt .. $], txmdmiddle, 1);

                    txStreamer.send(null, txmdend, 0.0001);
                    writeln("STOP");
                    continue;
                }
            }else if(!isSendMode){
                Thread.sleep(10.msecs);
            }
        }
    }
    catch(Throwable ex){
        writeln(ex);
    }
}


void receiveWorker(
    string addr,
    double rate,
    double gain,
    double freq,
    string clockSource,
    shared(Channel!ICA) channel
)
{
    try{
        auto usrp = USRP(MultiDeviceAddress([addr]));

        usrp.rxRate = rate;
        writeln("Actual RX Rate: ", usrp.rxRate);

        usrp.rxGain = gain;
        writeln("Actual RX Gain: ", usrp.rxGain);

        usrp.rxFreq =freq;
        writefln("Actual RX freq: %s [MHz]", usrp.rxFreq / 1e6);

        usrp.clockSource = clockSource;
        writefln("Actual clock source: %s", usrp.clockSource);

        // usrp.timeSource = clockSource;
        // writefln("Actual clock source: %s", usrp.timeSource);

        auto rxStreamer = usrp.makeRxStreamer(StreamArgs("fc32", "sc16", "", [0]));
        auto maxlen = rxStreamer.maxNumSamps;
        writefln("Buffer size in samples: %s", maxlen);

        rxStreamer.issue(StreamCommand.startContinuous());

        auto md = makeRxMetaData();

        auto buffer = new Complex!float[maxlen];

        bool lastPut = false;
        real avgpower = 0;
        size_t rxcnt;
        while(!endFlag)
        {
            auto nrecv = rxStreamer.recv(buffer, md, 0.1);
            if(nrecv == 0) continue;
            ++rxcnt;

            if(rxcnt > 1000 && rxcnt < 1100){
                avgpower += buffer[0 .. nrecv].map!"a.re^^2+a.im^^2".sum();
            }else if(rxcnt == 1100){
                avgpower /= 100;
            }
            auto pw = buffer[0 .. nrecv].map!"a.re^^2+a.im^^2".sum();
            auto snr = pw / avgpower;
            if(rxcnt % 1000 == 0 && snr > 2) writefln("Recv: %s", 10*log10(snr));
            if(md.errorCode != RxMetaData.ErrorCode.NONE)
            {
                writeln(md.errorCode);
                break;
            }

            if(snr > 2){
                if(!lastPut){
                    GC.disable();
                    channel.put(cast(ICA)null);
                }
                lastPut = true;
                channel.put(cast(ICA)(buffer[0 .. nrecv].dup));
            }else{
                if(lastPut){
                    channel.put(cast(ICA)null);
                    GC.enable();
                }
                lastPut = false;
            }
        }
    }
    catch(Throwable ex){
        writeln(ex);
    }
}


void txMain(Mod)(
    shared(Channel!ICA) txchannel,
    in size_t[] data,
    Mod mod,
    real payloadGain = 0.1,
    bool isSimulation = true,
)
{
    auto noiseGen = boxMullerNoise();

    OFDM!C ofdm = new OFDM!C(Constant.OFDM.nFFT, Constant.OFDM.nCP, Constant.OFDM.nSC, Constant.nOverSampling);

    ICA txsignal = ()
    {
        auto signal = modulateTX(ofdm, mod, cast(const(ubyte)[])data, payloadGain).dup;
        return signal.assumeUnique;
    }();


    writeln("POWER: ", txsignal[$/2 .. $].map!sqAbs.sum / (txsignal.length / 2));

    if(!isSimulation){
        Thread.sleep(5.seconds);
    }

    while(!endFlag){
        writefln("Send %s [bits]", data.length * typeof(data[0]).sizeof * 8);
        txchannel.put!ICA(null);
        txchannel.put!ICA(txsignal);
        txchannel.put!ICA(null);
        if(!isSimulation)
            Thread.sleep(500.msecs);
        else{
            Thread.sleep(50.msecs);
        }
    }
}


void rxMain(Mod)(
    shared(Channel!ICA) rxchannel,
    in size_t[] data,
    Mod mod,
    string paygain_dB_str,
    bool isSimulation = true,
)
{
    try{
        OFDM!C ofdm = new OFDM!C(Constant.OFDM.nFFT, Constant.OFDM.nCP, Constant.OFDM.nSC, Constant.nOverSampling);
        PreambleDetector detector = new PreambleDetector(preambles);

        File fileTx = File(format("transmitOFDM_%sdB_%s.dat", paygain_dB_str, Mod.stringof[0..4]), "w");
        File file = File(format("receivedOFDM_%sdB_%s.dat", paygain_dB_str, Mod.stringof[0..4]), "w");

        size_t recvCount;
        size_t errorCount;
        //size_t[2][] list;

        static struct Result
        {
            size_t err;
            size_t bits;
            ICA payload;
            ICA subcarriers;
        }

        Result[] list;

        while(!endFlag){
            if(auto p = rxchannel.pop!ICA)
            {
                //writefln("GET FROM RXWK: %s [size]", (*p).length);
                ICA received = *p;
                while(1){
                    if(auto q = rxchannel.pop!ICA){
                        //writefln("GET FROM RXWK: %s [size]", (*q).length);
                        if((*q).length){
                            received ~= *q;
                        }else{
                            break;
                        }
                    }
                }

                //writefln("Process: %s", received.length);

                received = findPreamble(detector, received);


                C[] receivedSubcarriers;
                ubyte[] bins = demodulateRX(ofdm, mod, received, receivedSubcarriers);

                immutable totalSC = data.length * typeof(data[0]).sizeof * 8 / Mod.symInputLength;
                immutable totalHeaderSamples = numOfHeaderSamples(ofdm, mod);
                immutable payloadReplica = modulateImpl(ofdm, mod, cast(const(ubyte[]))data).dup;
                immutable totalPayloadSamples = payloadReplica.length;

                if(bins.length && bins.length % size_t.sizeof == 0){
                    //assert(bins.length % size_t.sizeof == 0);
                    auto rcv = cast(const(size_t)[])bins;
                    if(rcv.length == data.length && receivedSubcarriers.length >= totalSC){
                        writefln("received!: %s [bits]", bins.length * 8);

                        immutable bits = bins.length * 8;

                        BitArray rxbits = BitArray(cast(void[])rcv, bits);
                        BitArray txbits = BitArray(cast(void[])data, bits);
                        size_t partErrorCount, partRecvCount;
                        foreach(i; 0 .. bits/2){    // 最初半分だけでテスト
                            if(rxbits[i] != txbits[i]){
                                ++errorCount;
                                ++partErrorCount;
                            }
                            ++recvCount;
                            ++partRecvCount;
                        }

                        writefln("%s: %s / %s = %s",
                            list.length,
                            errorCount, recvCount, errorCount * 1.0 / recvCount);

                        auto payloadsignal = received[totalHeaderSamples .. $][0 .. totalPayloadSamples];
                        list ~= Result(partErrorCount, partRecvCount, payloadsignal, receivedSubcarriers[0 .. data.length * typeof(data[0]).sizeof * 8 / Mod.symInputLength].dup);

                        if(recvCount >= (totalBits*2)){
                            endFlag = true;

                            // 外れ値を抜く
                            Result[] middle;
                            {
                                list.sort!"a.err*1.0 / a.bits < b.err * 1.0 / b.bits"();
                                middle = list[$/4 .. $ - $/4];
                                enforce(middle.map!"a.bits".sum() == totalBits);
                            }

                            auto sumerror = middle.map!"a.err".sum();
                            auto sumtotal = middle.map!"a.bits".sum();
                            writefln("BER: %s", sumerror * 1.0 / sumtotal);

                            //// ファイルに書き込む
                            //{
                            //    foreach(m; middle)
                            //        file.rawWrite(m.subcarriers);
                            //}
                            //// CSVファイルに先頭1000サンプル書き込む
                            //{
                            //    File csv = File(format("baseband_%s.csv", Mod.stringof[0..4]), "w");
                            //    foreach(e; middle.map!"a.subcarriers".joiner.take(1000))
                            //        csv.writefln("%s,%s", e.re, e.im);
                            //}

                            // ファイルにOFDM信号を書き出す
                            foreach(m; middle)
                                file.rawWrite(m.payload[$/2 .. $]);

                            fileTx.rawWrite(payloadReplica);

                            return;
                        }
                    }
                }
            }
            else
                Thread.sleep(100.msecs);
        }
    }catch(Throwable ex){
        writeln(ex);
        throw ex;
    }
}
