FROM --platform=linux/amd64 ubuntu:20.04 as builder

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential cmake clang

COPY . /repo
WORKDIR /repo/metis-5.1.0
RUN make config
RUN make -j8
RUN make install
RUN clang++ /repo/mayhem/fuzz_libmetis__InitRandom.cpp -fsanitize=fuzzer,address /repo/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a -o /fuzz

FROM ubuntu:20.04 as package
COPY --from=builder /usr/local/bin/graphchk /
COPY --from=builder /fuzz /fuzz_libmetis__InitRandom
