helloworld.exe: hello_world.o hello_world_cuda.o
	CC -o $@ $^

hello_world.o: hello_world.cpp hello_world.hpp
	CC -c -o $@ $<

hello_world_cuda.o: hello_world.cu hello_world.hpp
	nvcc -c -o $@ $<
