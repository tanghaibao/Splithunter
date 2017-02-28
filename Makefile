PROG=Splithunter
SEQ=SeqLib
CXX=g++
CXXFLAGS+=-I$(SEQ) -std=c++11 -DNDEBUG -W -Wall -pedantic -fopenmp
SEQLIB=$(SEQ)/bin/libseqlib.a
LDFLAGS+=$(SEQ)/bin/libseqlib.a $(SEQ)/bin/libbwa.a $(SEQ)/bin/libfml.a $(SEQ)/bin/libhts.a

OS := $(shell uname)
ifeq ($(OS), Darwin)
	LDFLAGS+=-lpthread -lz -lcurl -lcrypto
else
	LDFLAGS+=-lpthread -lrt -lz -lcurl -lcrypto
endif

SRCS=Splithunter.cpp
OBJS=$(SRCS:.cpp=.o)

default: all
all: $(PROG)

$(SEQLIB):
	git submodule update --init --recursive
	@cd SeqLib && ./configure LDFLAGS='-lcurl -lcrypto' && make -j 32 CXXFLAGS='-std=c++11' && make install

$(PROG): $(OBJS) $(SEQLIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cc
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

clean:
	@rm -f $(PROG) *.o

.PHONY: default all clean
