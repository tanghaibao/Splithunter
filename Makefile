PROG=Splithunter
SEQ=SeqLib
CXX=g++
CXXFLAGS+=-I$(SEQ) -I$(SEQ)/htslib -std=c++11 -DNDEBUG -W -Wall -pedantic
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
	./install.sh

$(PROG): $(OBJS) $(SEQLIB)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) -O3 -c -o $@ $^ $(CXXFLAGS)

clean:
	@rm -rf $(SEQ)/bin
	@rm -f $(PROG) *.o

.PHONY: default all clean
