PROGS=fpu

CPP=g++

# assembler options
ASM=-S -fverbose-asm  
#ASM+=-masm=intel -DASM_INTEL		# nicer asm intel syntax

# compile options
CFLAGS=-Wall -O3 -funroll-loops

# architechture (e.g. x86_32 or x86_64 bit, usage of sse2 etc)
ARCH=-march=native
#ARCH=-march=i386 -m32
#ARCH+=-mfpmath=sse -msse2 

#LFLAGS=-static

all:	${PROGS}


# cancelling default rules as they would short cut assembly and object stages
%.o : %.cpp

% : %.cpp

% : %.s


# cpp to assembly
%.s : %.cpp
	${CPP} ${ASM} ${CFLAGS} ${ARCH} $<

# assembly to object files
%.o : %.s
	${CPP} ${ARCH} -c $< 

# object to executable
%: %.o
	${CPP} ${ARCH} $< ${LFLAGS} -o $@

# mark .o and .s files as precious,
# otherwise would be deleted as considered intermediates
.PRECIOUS: %.o %.s


clean:
	rm -f *.s *.o ${PROGS} a.out
