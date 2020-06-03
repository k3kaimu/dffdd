# dffdd: Digital Filter library for Full-Duplex written by D

This library ``dffdd'' is a software library written by D programming language for baseband signal simulation of in-band full-duplex transceivers.
The author created this library to take the results of full-duplex studies[2,3].

This library provides following primitive modules and self-interference cancellers:

+ Phase Shift Keying modulator, demodulator
+ OFDM modulator, demodulator
+ Low Density Parity Check (LDPC) code
+ FIR and IIR filter
+ Memoryless Nonlinearity: Saleh, Rapp, etc...
+ Linear Least Squares methods
+ Adaptive Algorithms: LMS, NLMS, RLS
+ Parallel Hammerstein cancellers on time-domain[1] and frequency-domain[2]
+ Iterative Nonlinear canceller[3]

## How to use

In this instruction, we use [Homebrew](https://brew.sh/) to install compilers and libraries such as FFTW, OpenBLAS, and NLopt.
If you are using Linux, you can use the package manager of your distribution instead of Homebrew.
Of course, you can build an environment with Homebrew on Linux.

### Step 0: Install D compiler and tools

```sh
$ brew install dmd dub ldc
```

### Step 1: Install FFTW, OpenBLAS, and NLopt
This library requires `FFTW`, `OpenBLAS`, and `NLopt`.

```sh
$ brew install fftw openblas nlopt
```

### Step 2: Clone repositories

```sh
$ git clone https://github.com/k3kaimu/dffdd.git
$ git clone https://github.com/k3kaimu/carbon.git
```


### Step 3: Build and run a full-duplex simulator

~~~~
$ cd dffdd/apps/sim
$ dub --build=release
~~~~



## References:

+ [1]: D. Korpi, T. Huusari, Y. Choi, L. Anttila, S. Talwar and M. Valkama, "Digital self-interference cancellation under nonideal RF components: Advanced algorithms and measured performance," 2015 IEEE 16th International Workshop on Signal Processing Advances in Wireless Communications (SPAWC), Stockholm, 2015, pp. 286-290, doi: 10.1109/SPAWC.2015.7227045. 
+ [2]: K. Komatsu, Y. Miyaji and H. Uehara, "Basis Function Selection of Frequency-Domain Hammerstein Self-Interference Canceller for In-Band Full-Duplex Wireless Communications," in IEEE Transactions on Wireless Communications, vol. 17, no. 6, pp. 3768-3780, June 2018, doi: 10.1109/TWC.2018.2816061.
+ [3]: K. Komatsu, Y. Miyaji and H. Uehara, "Iterative Nonlinear Self-Interference Cancellation for In-Band Full-Duplex Wireless Communications Under Mixer Imbalance and Amplifier Nonlinearity," in IEEE Transactions on Wireless Communications, doi: 10.1109/TWC.2020.2983407.
