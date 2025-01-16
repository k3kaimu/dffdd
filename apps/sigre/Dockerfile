FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]

RUN <<EOF
    apt update
    apt install -y curl wget git xz-utils libxml2 gcc-11 liblapacke-dev libopenblas-dev libfftw3-dev
EOF

# install ldc
RUN <<EOF                                    
    cd /
    mkdir -p /dlang && wget https://dlang.org/install.sh -O /dlang/install.sh
    chmod +x /dlang/install.sh
    /dlang/install.sh install -p /dlang ldc-1.39.0
EOF

RUN <<EOF
    cd /root
    git clone https://github.com/k3kaimu/carbon.git
    git clone https://github.com/k3kaimu/lubeck.git
EOF

COPY . /root/build

RUN <<EOF
    cd /root/build/apps/sigre
    source /dlang/ldc-1.39.0/activate
    dub build --build=release --single sigre.d
EOF

RUN mkdir /bundle
RUN cp -r --parents \
    /root/build/apps/sigre/sigre \
    /lib/x86_64-linux-gnu/libopenblas.so* \
    /lib/x86_64-linux-gnu/liblapacke.so* \
    /lib/x86_64-linux-gnu/libfftw3.so* \
    /lib/x86_64-linux-gnu/libfftw3f.so* \
    /lib/x86_64-linux-gnu/libfftw3l.so* \
    /lib/x86_64-linux-gnu/libm.so* \
    /lib/x86_64-linux-gnu/libgcc_s.so* \
    /lib/x86_64-linux-gnu/libc.so* \
    /lib/x86_64-linux-gnu/libgfortran.so* \
    /lib/x86_64-linux-gnu/libblas.so* \
    /lib/x86_64-linux-gnu/liblapack.so* \
    /lib/x86_64-linux-gnu/libtmglib.so* \
    /bundle

RUN mkdir /alternatives
RUN cp \
    /etc/alternatives/libblas.so* \
    /etc/alternatives/liblapack.so* \
    /etc/alternatives/libopenblas.so* \
    /alternatives

FROM ubuntu:24.04
COPY --from=0 /bundle /bundle
COPY --from=0 /alternatives /etc/alternatives
ENV LD_LIBRARY_PATH=/bundle/lib/x86_64-linux-gnu
WORKDIR /work
ENTRYPOINT ["/bundle/root/build/apps/sigre/sigre"]