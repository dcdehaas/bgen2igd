FROM ubuntu:22.04

RUN apt update && \
    apt install -y git build-essential cmake

COPY . /bgen2igd

RUN cd /bgen2igd && mkdir cpp_build && cd cpp_build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j && \
    cp bgen2igd /usr/bin/
