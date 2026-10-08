all: usp-jewel usp-jewel-inter

# path to LHAPDF library, write yours here
LHAPDF_PATH := /home/leonardo/Applications/LHAPDF/lib 

FC := gfortran
FFLAGS := -O2

usp-jewel: jewel260mod-usp.o medium-hydro.o reader-hydro.o interpolate.o subtraction.o pythia6425mod-lhapdf6.o meix.o
	$(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF -lstdc++

clean:
	rm -f jewel*.o
	rm -f medium-*.o interpolate.o reader-*.o 
	rm -f pythia6425mod-lhapdf6.o meix.o subtraction.o
	rm -f *~

.PHONY: all
