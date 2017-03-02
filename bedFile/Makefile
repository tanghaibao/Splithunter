CXX = g++
OBJ_DIR = ./
# -------------------
# define our includes
# -------------------
INCLUDES = -I.

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= bedFile.cpp gzstream.cpp fileType.cpp
OBJECTS= $(SOURCES:.cpp=.o)
BEDLIB = libbedFile.a

$(BEDLIB): $(OBJECTS)
	ar -rcs $@ $<

%.o: %.cpp
	@echo "  * Compiling" $(*F).cpp
	$(CXX) -O3 -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES)

default: all
all: $(OBJECTS)

clean:
	@echo "Cleaning up."
	@rm -f *.o *.a

.PHONY: clean
