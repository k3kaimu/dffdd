# dffdd: Digital Filter library for Full-Duplex written by D

This library ``dffdd'' is a software library written by D programming language for baseband signal simulation of in-band full-duplex transceivers.
This library provides following primitive modules and self-interference cancellers:

+ Phase Shift Keying modulator, demodulator
+ OFDM modulator, demodulator
+ Low Density Parity Check (LDPC) code
+ FIR and IIR filter
+ Memoryless Nonlinearity: Saleh, Rapp, etc...
+ Linear Least Squares methods
+ Adaptive Algorithms: LMS, NLMS, RLS
+ Parallel Hammerstein cancellers on time-domain and frequency-domain
+ Iterative Nonlinear canceller

## How to use

In this instruction, we use [Homebrew](https://brew.sh/) to install compilers, and libraries such as FFTW, OpenBLAS, NLopt.
If you are using Linux, you can use the package manager of your distribution instead of Homebrew.
Of course, you can build an environment with Homebrew on Linux.

### Step 0: Install D compiler and tools

```sh
$ brew install dmd dub
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


## Step 3: Build and run a full-duplex simulator

~~~~
$ cd apps/sim
$ dub --build=release
~~~~
