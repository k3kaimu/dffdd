FROM alpine

RUN apk add --no-cache bash curl wget git build-base lapack-dev openblas-dev fftw-dev ldc dub gcc

RUN <<EOF
    cd /root
    git clone https://github.com/k3kaimu/carbon.git
    git clone https://github.com/k3kaimu/lubeck.git
EOF

COPY . /root/build

RUN <<EOF
    cd /usr/bin/
    ln -s gcc gcc-11
    cd /root/build/apps/sigre
    dub build --build=release --single sigre.d
    cp ./sigre /usr/bin
EOF

RUN mkdir /bundle
RUN cp -r --parents \
    /usr/bin/sigre \
    /usr/lib/libopenblas.so* \
    /usr/lib/libopenblasp-r*.so* \
    /usr/lib/libfftw3.so* \
    /usr/lib/libfftw3f.so* \
    /usr/lib/libfftw3l.so* \
    /usr/lib/liblapacke.so* \
    /usr/lib/libphobos2-ldc-shared.so* \
    /usr/lib/libdruntime-ldc-shared.so* \
    /usr/lib/libgcc_s.so* \
    /usr/lib/libgfortran.so* \
    /usr/lib/libunwind.so* \
    /usr/lib/libquadmath.so* \
    /bundle

FROM alpine
COPY --from=0 /bundle /
WORKDIR /work
ENTRYPOINT ["/usr/bin/sigre"]