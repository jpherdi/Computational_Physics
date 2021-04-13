CXX       := g++
CXX_FLAGS := -std=c++17 -ggdb

BIN     := bin
SRC     := src
INCLUDE := include
EIGEN := ~/Downloads/eigen/

LIBRARIES   := #lib
EXECUTABLE  := main


all: $(BIN)/$(EXECUTABLE)

run: clean all
	clear
	./$(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	$(CXX) $(CXX_FLAGS) -I $(INCLUDE) -I $(EIGEN) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*