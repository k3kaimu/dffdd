module spec;


struct Constant
{
    enum uint nOverSampling = 4;
    enum real samplingFreq = 4e6;
    enum bool isFDMode = false;
    enum real centerFreq = 5.015e9;

    struct OFDM
    {
        enum size_t numOfBitsOfSC = 2;
        enum size_t nFFT = 64;
        enum size_t nCP = 16;
        enum size_t nSC = 52;
        enum size_t numOfTotalSamples = (nFFT + nCP) * nOverSampling;
    }

    struct PreambleDetector
    {
        enum size_t numOfPreamble = 1023 * nOverSampling;
        enum size_t numOfChunk = 1024 * nOverSampling;
        enum real acqTH = 10;
    }

    struct EST
    {
        enum size_t numOfESTCycle = 10; // 100Hzの誤差を許容する
        enum size_t numOfTotalSamples = PreambleDetector.numOfPreamble * numOfESTCycle;
    }

    struct RTSCTS
    {
        enum size_t numOfLearningSymbol = 12;
        enum size_t numOfTotalSamples = PreambleDetector.numOfPreamble + OFDM.numOfTotalSamples * numOfLearningSymbol;
    }

    struct DAT
    {
        enum size_t numOfPayloadSymbols = 128 + 7;
        enum size_t numOfPayloadSamples = OFDM.numOfTotalSamples * numOfPayloadSymbols;
        enum size_t numOfTotalSamples = PreambleDetector.numOfPreamble + numOfPayloadSamples;
    }

    struct ACK
    {
        enum size_t numOfTotalSamples = 1023 * nOverSampling;
    }
}
