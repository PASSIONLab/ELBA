K?=31
L?=15
U?=35
COMPILE_TIME_PARAMETERS=-DKMER_SIZE=$(K) -DLOWER_KMER_FREQ=$(L) -DUPPER_KMER_FREQ=$(U)

MPICH=/usr/local/Cellar/mpich/4.1.1
MPICH_INC=-I$(MPICH)/include
MPICH_LIB=-L$(MPICH)/lib
MPICH_FLAGS=
FLAGS=$(COMPILE_TIME_PARAMETERS) -O2 -Wno-maybe-uninitialized -Wno-deprecated -std=c++17 -I./include -I./src

COMBBLAS=./CombBLAS
COMBBLAS_INC=$(COMBBLAS)/include/CombBLAS
COMBBLAS_SRC=$(COMBBLAS)/src
INCADD=-I$(COMBBLAS)/include/ -I$(COMBBLAS)/psort-1.0/include/ -I$(COMBBLAS)/usort/include/ -I$(COMBBLAS)/graph500-1.2/generator/include/

UNAME_S:=$(shell uname -s)

ifeq ($(UNAME_S),Linux)
COMPILER=CC
else ifeq ($(UNAME_S),Darwin)
COMPILER=g++-12
FLAGS+=$(MPICH_INC)
MPICH_FLAGS+=$(MPICH_LIB) -L/usr/local/opt/libevent/lib -lmpi
endif

OBJECTS=obj/Logger.o \
		obj/FastaIndex.o \
		obj/DistributedFastaData.o \
		obj/DnaSeq.o \
		obj/DnaBuffer.o \
		obj/HashFuncs.o \
		obj/HyperLogLog.o \
		obj/Bloom.o \
		obj/KmerOps.o \
		obj/ReadOverlap.o \
		obj/CommGrid.o \
		obj/MPIType.o

all: elba

test: elba
	./runtests.sh

elba: obj/main.o $(OBJECTS)
	@echo CXX -c -o $@ $^
	@$(COMPILER) $(FLAGS) $(INCADD) -o $@ $^ $(MPICH_FLAGS) -lz

obj/%.o: src/%.cpp
	@echo CXX $(COMPILE_TIME_PARAMETERS) -c -o $@ $<
	@$(COMPILER) $(FLAGS) $(INCADD) -c -o $@ $<

obj/main.o: src/main.cpp include/common.h src/Kmer.cpp include/Kmer.hpp src/KmerOps.cpp include/KmerOps.hpp src/ReadOverlap.cpp include/ReadOverlap.hpp include/KmerIntersect.hpp
obj/Logger.o: src/Logger.cpp include/Logger.hpp
obj/FastaIndex.o: src/FastaIndex.cpp include/FastaIndex.hpp
obj/DistributedFastaData.o: src/DistributedFastaData.cpp include/DistributedFastaData.hpp
obj/KmerOps.o: src/KmerOps.cpp include/KmerOps.hpp
obj/ReadOverlap.o: src/ReadOverlap.cpp include/ReadOverlap.hpp
obj/DnaSeq.o: src/DnaSeq.cpp include/DnaSeq.hpp
obj/DnaBuffer.o: src/DnaBuffer.cpp include/DnaBuffer.hpp
obj/HashFuncs.o: src/HashFuncs.cpp include/HashFuncs.hpp

obj/CommGrid.o: $(COMBBLAS_SRC)/CommGrid.cpp $(COMBBLAS_INC)/CommGrid.h
	@echo CXX -c -o $@ $<
	@$(COMPILER) $(FLAGS) $(INCADD) -c -o $@ $<

obj/MPIType.o: $(COMBBLAS_SRC)/MPIType.cpp $(COMBBLAS_INC)/MPIType.h
	@echo CXX -c -o $@ $<
	@$(COMPILER) $(FLAGS) $(INCADD) -c -o $@ $<

clean:
	rm -rf *.o obj/*.o *.dSYM *.out *.mtx

gitclean: clean
	git clean -f
