#!/bin/bash

cd SeqLib
cd htslib && ./configure --enable-libcurl && cd ..
./configure LDFLAGS='-lcurl -lcrypto'
make CXXFLAGS='-std=c++11' -j 32
make install
