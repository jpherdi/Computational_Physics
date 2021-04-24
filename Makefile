CXXFLAGS = -O3 -ffast-math -march=native -std=c++17 

all: profilerExample.cpp profiler.cpp profiler.h
	g++ $(CXXFLAGS) profilerExample.cpp profiler.cpp -o profilerExample
	./profilerExample
	
# 'make clean' entfernt Zwischendateien
clean:
	rm -f ./profilerExample

# .PHONY gibt an, dass clean auch funktioniert, wenn eine Datei namens 'clean' existiert
.PHONY: clean