import std.stdio;
import std.stdio;
import std.getopt;
import std.range;
import std.algorithm;
import std.array;
import std.format;
import std.string;
import std.numeric;
import std.math;
import std.complex;
import std.datetime;

import carbon.stream;

enum size_t P = 8;
enum size_t N = 8;
enum size_t Total = 1024*64;
enum real beta = 0.00001;
enum real sampFreq = 400e3;
enum size_t blockSize = 1024;
enum real PLLB = 20;
enum real PLL_K = PLLB / 0.53;
enum real PLLAW = 1.414 * PLL_K;
enum real PLLW2 = PLL_K ^^ 2;
enum real FLLB = 250;
enum real FLLW = FLLB / 0.25;

private F atan(F)(F y, F x)
{
    return x.signbit ? atan2(-y, -x) : atan2(y, x);
}

void main(string[] args)
{
    string sendFilename, recvFilename, outFilename;
    bool bSpeedMode = false;

    getopt(args,
        "sendData", &sendFilename,
        "recvData", &recvFilename,
        "outData", &outFilename,
        "speedMode", &bSpeedMode
    );

    File sendFile = File(sendFilename),
         recvFile = File(recvFilename),
         outFile = File(outFilename, "w");

    recvFile.seek(919 * 8);

    creal[P][] w = new creal[P][](N),
               x = new creal[P][](N);

    foreach(ref e; w) e[] = 0+0i;
    foreach(ref e; x) e[] = 0+0i;

    cfloat[] sendBuf = new cfloat[blockSize],
              recvBuf = new cfloat[blockSize],
              writeBuf = new cfloat[blockSize];

    double[] fftResultRecv = new double[blockSize],
             fftResultSIC = new double[blockSize];

    //creal oldPrompt = 0+0i;
    //real oldCarrErr = 0;
    //real carrNCOFreq = 400;
    //immutable DT = blockSize / sampFreq;
    real pmaxRecv = 1E-50L, pmaxSend = 1E-300L;

    Fft fftObj = new Fft(blockSize);

    //auto nco = lutNCO!(std.math.expi, 1024*16)(carrNCOFreq, DT, 0);
    auto startTime = Clock.currTime;

    foreach(_unused_; 0 .. Total)
    {
        auto sendGets = sendFile.rawRead(sendBuf);
        auto recvGets = recvFile.rawRead(recvBuf);
        assert(sendGets.length == sendBuf.length);
        assert(recvGets.length == recvBuf.length);

        if(_unused_ == 0) foreach(_, n; sendGets)
        {
            auto pSend = n.abs()^^2,
                 pRecv = recvGets[_].abs()^^2;

            if(pRecv > pmaxRecv) pmaxRecv = pRecv;
            if(pSend > pmaxSend) pmaxSend = pSend;
        }

        creal prompt = 0+0i;
        real  pmaxRecvMax = 0, pmaxSendMax = 0;

        foreach(_, ref n; sendGets){
            immutable pSend = n.abs()^^2,
                      pRecv = recvGets[_].abs()^^2;

            if(pRecv > pmaxRecvMax) pmaxRecvMax = pRecv;
            if(pSend > pmaxSendMax) pmaxSendMax = pSend;

            //recvGets[_] *= nco.front;
            //nco.popFront();

            /* for PLL */
            n *= sqrt(1 / pmaxSend);
            recvGets[_] *= sqrt(1 / pmaxRecv);
            prompt += n.conj * recvGets[_];

            immutable nabs = n.abs;
            foreach(p, ref e; x[0]){
                e = (p+1)%2==0 ? 0+0i : 1+0i;
                foreach(__; 0 .. p)
                    e *= nabs;
                e *= n;
            }

            creal d = recvGets[_];
            creal y = 0+0i;
            foreach(i; 0 .. N){
                foreach(p; 0 .. P)
                    y += w[i][p] * x[i][p];
            }

            creal e = d - y;
            foreach_reverse(i; 0 .. N){
                foreach(p; 0 .. P)
                    w[i][p] += beta * e * x[i][p].conj;

                if(i != 0) x[i] = x[i-1];
            }

            if(!bSpeedMode && _unused_ % 1000 == 0)
            {
                outFile.writefln("%s,%s,%s,%s,%s,%s,%s,%s,%s,", _unused_, n.re, n.im, recvGets[_].re, recvGets[_].im, y.re, y.im, e.re, e.im);
            }

            writeBuf[_] = e;
        }
        outFile.flush();

        pmaxRecv = (pmaxRecv + pmaxRecvMax) / 2;
        pmaxSend = (pmaxSend + pmaxSendMax) / 2;

        ///* update PLL */
        //immutable carrErr = atan(prompt.im, prompt.re) / 2 / PI,
        //          freqErr = atan2(oldPrompt.re * prompt.im - prompt.re * oldPrompt.im,
        //                          abs(oldPrompt.re * prompt.re) + abs(oldPrompt.im * prompt.im)) / PI;

        //carrNCOFreq += PLLAW * (carrErr - oldCarrErr)
        //             + PLLW2 * DT * carrErr
        //             + FLLW * DT * freqErr;

        //nco.freq = carrNCOFreq;

        //oldPrompt = prompt;
        //oldCarrErr = carrErr;

        // fft
        if(!bSpeedMode){
            auto cpxBuf1 = new Complex!double[blockSize];
            auto cpxBuf2 = writeBuf.map!(a => complex!double(a.re, a.im)).array();
            fftObj.fft(cpxBuf2, cpxBuf1);
            cpxBuf2[0 .. $/2] = cpxBuf1[$/2 .. $];
            cpxBuf2[$/2 .. $] = cpxBuf1[0 .. $/2];

            auto recvBuf1 = new Complex!double[blockSize];
            auto recvBuf2 = recvBuf.map!(a => complex!double(a.re, a.im)).array();
            fftObj.fft(recvBuf2, recvBuf1);
            recvBuf2[0 .. $/2] = recvBuf1[$/2 .. $];
            recvBuf2[$/2 .. $] = recvBuf1[0 .. $/2];

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] += recvBuf2[i].re^^2 + recvBuf2[i].im^^2 / 1000;
                fftResultSIC[i] = cpxBuf2[i].re^^2 + cpxBuf2[i].im^^2 / 1000;
            }
        }

        if(!bSpeedMode && _unused_ % 1000 == 0){
            foreach(i; 0 .. blockSize)
            {
                auto before = 10*log10(fftResultRecv[i]),
                     after = 10*log10(fftResultSIC[i]);

                outFile.writefln("%s,%s,%s,%s,%s,%s,", _unused_, i, i*400e3/blockSize-200.0e3, before, after, before - after);
            }
            outFile.flush();

            foreach(i; 0 .. blockSize){
                fftResultRecv[i] = 0;
                fftResultSIC[i] = 0;
            }
        }

        writefln("%s     %s       %s[k samples/s]", _unused_, writeBuf[$-1].abs, (_unused_+1) * blockSize / ((Clock.currTime - startTime).total!"msecs"() / 1000.0L) / 1000.0L);
    }
}
