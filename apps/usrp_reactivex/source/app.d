module app;

import std.stdio;

import uhd.usrp;
import uhd.capi.types.metadata;
import carbon.channel;
import std.random;
import std.complex;
import std.concurrency;
import std.algorithm;
import std.range;
import core.thread;
import std.math;
import std.datetime;
import std.bitmanip;

import dffdd.mod.ofdm;
import dffdd.mod.qpsk;
import dffdd.dsp.convolution;
import dffdd.utils.fft;

import zmqd;
import preamble;
import spec;
import impl;

import core.memory;

shared bool endFlag;



void main(string[] args)
{
    auto usrpAddr = args[1];
    auto clockSource = args[2];
    auto txport = args[3];
    auto rxport = args[4];

    // writeln("Edit source/app.d to start your project.");
    auto txChannel = channel!ICA();
    auto rxChannel = channel!ICA();

    // if(args[3] == "tx")
    auto txtid = spawn(
                    &transmitWorker,
                    usrpAddr,
                    1e6,
                    20,
                    5.015e9,
                    clockSource,
                    txChannel);
    
    // if(args[3] == "rx")
    auto rxtid = spawn(
                    &receiveWorker,
                    usrpAddr,
                    1e6,
                    35,
                    5.015e9,
                    clockSource,
                    rxChannel);

    auto zmqTX = spawn(
                    &zmqTXWorker,
                    txport,
                    txChannel);

    auto zmqRX = spawn(
                    &zmqRXWorker,
                    rxport,
                    rxChannel);

    auto txSender = Socket(SocketType.push);
    txSender.bind("tcp://*:" ~ txport);

    foreach(line; stdin.byLine){
        if(line.startsWith("s")){
            txSender.send("Hello, OKOKOKOKKKOKOKOKOKO");
        }else if(line.startsWith("e")){
            endFlag = true;
            break;
        }
    }
}

/+
void transmitWorker(
    string addr,
    double rate,
    double gain,
    double freq,
    string clockSource,
    shared(Channel!ICA) channel)
{
    auto usrp = USRP(MultiDeviceAddress([addr]));

    usrp.txRate = rate;
    usrp.rxRate = rate;
    writeln("Actual TX Rate: ", usrp.txRate);
    writeln("Actual RX Rate: ", usrp.rxRate);

    usrp.txGain = gain;
    usrp.rxGain = gain;
    writeln("Actual TX Gain: ", usrp.txGain);
    writeln("Actual RX Gain: ", usrp.rxGain);

    usrp.txFreq = freq;
    usrp.rxFreq = freq;
    writefln("Actual TX freq: %s [MHz]", usrp.txFreq / 1e6);
    writefln("Actual RX freq: %s [MHz]", usrp.rxFreq / 1e6);

    usrp.clockSource = clockSource;
    writefln("Actual clock source: %s", usrp.clockSource);

    auto txStreamer = usrp.makeTxStreamer(StreamArgs("fc32", "sc16", "", [0]));
    auto txMaxlen = 1024*1024;
    writefln("TXBuffer size in samples: %s", txMaxlen);

    auto rxStreamer = usrp.makeRxStreamer(StreamArgs("fc32", "sc16", "", [0]));
    auto rxMaxlen = rxStreamer.maxNumSamps;
    writefln("RXBuffer size in samples: %s", rxMaxlen);

    rxStreamer.issue(StreamCommand.startContinuous());

    // auto txmd = TxMetaData(false, 0, 0.1, true, false);
    auto rxmd = makeRxMetaData();

    auto txmdstart = TxMetaData(false, 0, 0.1, true, false);
    auto txmdmiddle = TxMetaData(false, 0, 0.1, false, false);
    auto txmdend = TxMetaData(false, 0, 0.1, false, true);


    auto buffer = new Complex!float[rxMaxlen];
    // foreach(i, ref e; )

    auto startTime = Clock.currTime;
    size_t rxelems;
    size_t rxcnt;
    real avgpower = 0;
    while(!endFlag){
        if(auto p = channel.pop!ICA){
            rxStreamer.issue(StreamCommand.stopContinuous());

            StopWatch sw;
            sw.start();
            txStreamer.send(null, txmdstart, 0.1);
            foreach(chunk; [(*p)]){
                size_t cnt;
                while(cnt < chunk.length)
                    cnt += txStreamer.send(chunk[cnt .. $], txmdmiddle, 0.0001);
            }
            txStreamer.send(null, txmdend, 0.1);
            sw.stop();
            writefln("TX Rate: %s [Msps]", (*p).length / (sw.peek.usecs * 1.0L / 1E6) );
            rxStreamer.issue(StreamCommand.startContinuous());
            startTime = Clock.currTime;
            rxelems = 0;
        }else{
            ++rxcnt;
            auto nrecv = rxStreamer.recv(buffer, rxmd, 0.1);
            if(rxcnt > 1000 && rxcnt < 1100){
                avgpower += buffer[0 .. nrecv].map!"a.re^^2+a.im^^2".sum();
            }else if(rxcnt == 1100){
                avgpower /= 100;
            }
            rxelems += nrecv;
            auto pw = buffer[0 .. nrecv].map!"a.re^^2+a.im^^2".sum();
            auto snr = 10*log10(pw / avgpower);
            if(rxcnt % 100 == 0 && snr > 3){
                auto nowTime = Clock.currTime;
                writefln("SNR: %s [dB], %s [Msps]", snr, rxelems / ((nowTime - startTime).total!"usecs" / 1.0E6)  );
            }
            if(rxmd.errorCode != uhd_rx_metadata_error_code_t.UHD_RX_METADATA_ERROR_CODE_NONE)
            {
                writeln(rxmd.errorCode);
                break;
            }
        }
            // Thread.sleep(1.msecs);
    }
}
+/


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
        auto txMaxlen = max(min(txStreamer.maxNumSamps, 1024), 1024);
        writefln("TXBuffer size in samples: %s", txMaxlen);

        auto txmdstart = TxMetaData(false, 0, 0.1, true, false);
        auto txmdmiddle = TxMetaData(false, 0, 0.1, false, false);
        auto txmdend = TxMetaData(false, 0, 0.1, false, true);

        while(!endFlag){
            if(auto p = channel.pop!ICA){
                writeln("START");
                txStreamer.send(null, txmdstart, 0.1);
                foreach(chunk; [*p]){
                    size_t cnt;
                    while(cnt < chunk.length){
                        cnt += txStreamer.send(chunk[cnt .. $], txmdmiddle, 0.1);
                        // writefln("%s : %s", cnt, chunk.length);
                    }
                }
                txStreamer.send(null, txmdend, 0.1);
                writeln("END");
            }else{
                Thread.sleep(100.msecs);
            }
                // Thread.sleep(1.msecs);
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
            ++rxcnt;
            auto nrecv = rxStreamer.recv(buffer, md, 0.1);
            if(rxcnt > 1000 && rxcnt < 1100){
                avgpower += buffer[0 .. nrecv].map!"a.re^^2+a.im^^2".sum();
            }else if(rxcnt == 1100){
                avgpower /= 100;
            }
            auto pw = buffer[0 .. nrecv].map!"a.re^^2+a.im^^2".sum();
            auto snr = pw / avgpower;
            if(rxcnt % 100 == 0 && snr > 2) writefln("Recv: %s", 10*log10(snr));
            if(md.errorCode != uhd_rx_metadata_error_code_t.UHD_RX_METADATA_ERROR_CODE_NONE)
            {
                writeln(md.errorCode);
                break;
            }

            if(snr > 2){
                if(!lastPut) GC.disable();
                lastPut = true;
                channel.put(cast(ICA)(buffer[0 .. nrecv].dup));
            }else{
                if(lastPut) GC.enable();
                lastPut = false;
            }
        }
    }
    catch(Throwable ex){
        writeln(ex);
    }
}


void zmqTXWorker(
    string zmqTXPort,
    shared(Channel!ICA) txchannel
)
{
    OFDM!C ofdm = new OFDM!C(64, 16, 52, 4);
    QPSK qpsk;

    auto txbinary = Socket(SocketType.pull);
    txbinary.connect("tcp://127.0.0.1:" ~ zmqTXPort);

    while(!endFlag){
        Frame frame = Frame();
        if(txbinary.receive(frame) != 0){
            ICA signal = modulateTX(ofdm, qpsk, frame.data).dup;
            txchannel.put!ICA(signal);
        }
    }
}


void zmqRXWorker(
    string zmqRXPort,
    shared(Channel!ICA) rxchannel
)
{
    OFDM!C ofdm = new OFDM!C(64, 16, 52, 4);
    QPSK qpsk;
    PreambleDetector detector = new PreambleDetector(preambles);


    auto rxbinary = Socket(SocketType.push);
    rxbinary.bind("tcp://*:" ~ zmqRXPort);

    while(!endFlag){
        if(auto p = rxchannel.pop!ICA)
        {
            ICA received = *p;
            bool successLast = true;
            while(1){
                if(auto q = rxchannel.pop!ICA){
                    successLast = true;
                    received ~= *q;
                }else{
                    if(!successLast) break;
                    successLast = false;
                    Thread.sleep(1.msecs);
                }
            }

            received = findPreamble(detector, received);
            auto bins = demodulateRX(ofdm, qpsk, received);
            if(bins.length){
                writefln("received!: %s", bins.length);
                rxbinary.send(bins);
            }
        }
        else
            Thread.sleep(100.msecs);
    }
}
