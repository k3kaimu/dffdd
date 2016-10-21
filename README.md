# dffdd: Digital Filter library for Full-Duplex written by D

## How to use

At first, you must install dmd and dub.

~~~~~
$ brew install dmd
$ brew install dub
~~~~~

This library requires `FFTW` and `OpenBLAS`.
They can be install by `homebrew`.

~~~~~
$ brew install fftw
$ brew install homebrew/science/openblas
~~~~~

In addition, you should clone this library from gitlab, `dranges` and `carbon` from github.

~~~~~
$ mkdir workdir && workdir
$ git clone https://whale0.comm.ee.tut.ac.jp/gitlab/komatsu/dffdd.git
$ git clone https://github.com/k3kaimu/carbon.git
$ git clone https://github.com/k3kaimu/dranges.git
~~~~~


## Example: build and run a full-duplex simulator

~~~~
$ cd apps/sim
$ dub --build=release
~~~~
