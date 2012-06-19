GCC = g++

simulation: simulation.cpp
	$(GCC) -o simulation simulation.cpp -O2
