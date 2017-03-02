SEQ=SeqLib
BED=bedFile
CXX=g++
CXXFLAGS+=-I. -I$(SEQ) -I$(SEQ)/htslib -std=c++11 -DNDEBUG -W -Wall -pedantic
SEQLIB=$(SEQ)/bin/libseqlib.a
LDFLAGS+=$(SEQ)/bin/libseqlib.a $(SEQ)/bin/libbwa.a $(SEQ)/bin/libfml.a $(SEQ)/bin/libhts.a $(BED)/libbedFile.a

OS := $(shell uname)
ifeq ($(OS), Darwin)
	LDFLAGS+=-lpthread -lz -lcurl -lcrypto
else
	LDFLAGS+=-lpthread -lrt -lz -lcurl -lcrypto
endif

SRCS=Splithunter.cpp
OBJS=$(SRCS:.cpp=.o)

default: all
all: Splithunter BuildDB

$(SEQLIB):
	git submodule update --init --recursive
	./install.sh

Splithunter: Splithunter.o $(SEQLIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

BuildDB: BuildDB.o $(SEQLIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) -O3 -c -o $@ $^ $(CXXFLAGS)

clean:
	@echo "Cleaning up."
	@rm -rf $(SEQ)/bin
	@rm -f $(PROG) *.o

.PHONY: default all clean
