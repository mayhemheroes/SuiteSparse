FROM --platform=linux/amd64 ubuntu:22.04 as builder

# Configure apt with retries for transient network failures
RUN echo "Acquire::Retries \"5\";" > /etc/apt/apt.conf.d/80retries && \
    apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential cmake clang git

# Download and build METIS
RUN git clone --depth=1 https://github.com/KarypisLab/GKlib.git /gklib && \
    git clone --depth=1 --branch v5.1.1-DistDGL-v0.5 https://github.com/KarypisLab/METIS.git /metis-src

RUN cd /metis-src && \
    mkdir -p build && cd build && \
    cmake .. -DGKLIB_PATH=/gklib -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_C_COMPILER=clang -DSHARED=0 && \
    make -j8

COPY . /repo
RUN clang++ /repo/mayhem/fuzz_libmetis__InitRandom.cpp -fsanitize=fuzzer,address \
    /metis-src/build/libmetis/libmetis.a \
    -I /metis-src/include \
    -I /gklib \
    -o /fuzz

FROM ubuntu:22.04 as package
COPY --from=builder /fuzz /fuzz_libmetis__InitRandom
