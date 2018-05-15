# dffdd: Digital Filter library for Full-Duplex written by D

## How to use

At first, you must install dmd and dub.

~~~~~
$ brew install dmd
$ brew install dub
~~~~~

This library requires `FFTW`, `OpenBLAS`, and `NLopt`.
They can be installed by `homebrew`.

~~~~~
$ brew install fftw openblas nlopt
~~~~~

In addition, you should clone this library and carbon from gitlab and github, respectively.

~~~~~
$ mkdir workdir && cd workdir
$ git clone https://whale0.comm.ee.tut.ac.jp/gitlab/komatsu/dffdd.git
$ git clone https://github.com/k3kaimu/carbon.git
$ git clone https://github.com/k3kaimu/dranges.git
~~~~~


## Example: build and run a full-duplex simulator

~~~~
$ cd apps/sim
$ dub --build=release
~~~~
