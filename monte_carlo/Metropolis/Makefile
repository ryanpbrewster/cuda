pi.exe: metropolis.o mcpi.o
	CC -o $@ $^
metropolis.o: metropolis.cpp metropolis.hpp
	CC -c -o $@ $<
mcpi.o: mcpi.cu metropolis.hpp
	nvcc -c -o $@ $<

.PHONY: clean
clean:
	rm -f pi.exe metropolis.o mcpi.o
