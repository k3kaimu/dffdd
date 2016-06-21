module constant;

import dffdd.blockdiagram.noise : noisePower;

struct ConstantList
{ static: immutable:
    struct ConstExample1
    { static: immutable:
        real samplingFreq = 20e6 * 4;


        struct QAM
        { static: immutable:
            uint arity = 16;
        }


        struct OFDM
        { static: immutable:
            uint numOfFFT = 64;
            uint numOfCP = 16;
            uint numOfSubcarrier = 52;
            uint scaleOfUpSampling = 4;
            real PAPR = 10;                 // 10dB
        }


        struct ThermalNoise
        { static: immutable:
            real temperature = 300;
            real power = noisePower(samplingFreq / UpDownSampler.scaleOfUpSampling / OFDM.numOfFFT * OFDM.numOfSubcarrier,
                                    ThermalNoise.temperature);
        }


        struct UpDownSampler
        { static: immutable:
            uint scaleOfUpSampling = 8;
            real[] decimationFIRFilterTaps = [
                0.0064217,0.0030971,0.003501,0.0036559,0.0034831,0.0029178,0.0019183,0.00047411,-0.0013877,-0.003598,-0.0060443,-0.0085724,
                -0.010993,-0.013088,-0.014628,-0.015381,-0.015132,-0.0137,-0.010952,-0.0068179,-0.0013004,0.0055175,0.013472,0.022324,
                0.031764,0.041433,0.050934,0.059859,0.06781,0.074421,0.079383,0.082458,0.0835,0.082458,0.079383,0.074421,0.06781,
                0.059859,0.050934,0.041433,0.031764,0.022324,0.013472,0.0055175,-0.0013004,-0.0068179,-0.010952,-0.0137,-0.015132,
                -0.015381,-0.014628,-0.013088,-0.010993,-0.0085724,-0.0060443,-0.003598,-0.0013877,0.00047411,0.0019183,0.0029178,
                0.0034831,0.0036559,0.003501,0.0030971,0.0064217];

            size_t numOfShift = 67;
            //size_t numOfShift = 0;
        }


        struct PowerSpectralDensityAnalyzer
        { static: immutable:
            uint resolution = 1024;
        }


        struct BERCounter
        { static: immutable:
            ulong totalBits = 10_000;
        }


        struct TXIQMixer
        { static: immutable:
            real IIR = 25;  // 25dB
        }


        struct RXIQMixer
        { static: immutable:
            real IIR = 25;  // 25dB
        }


        struct PA 
        { static: immutable:
            real GAIN = 30;     // 30dB
            real IIP3 = 10;     // 10dBm
            real TX_POWER = 30; // 20dBm
        }

        struct LNA
        { static: immutable:
            real GAIN = 20;     // 20dB
            real NF = 4;        // 4dB
        }


        struct Quantizer
        { static: immutable:
            uint numOfBits = 14;
        }


        struct AWGN
        { static: immutable:
            real SNR = 24;      // 24dB
        }
    }
}
