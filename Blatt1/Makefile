CXX       := g++
CXX_FLAGS := -std=c++17 -ggdb

BIN     := bin
SRC     := src
EIGEN := ~/Downloads/eigen/

LIBRARIES   := #lib
EXECUTABLE  := main


all: $(BIN)/$(EXECUTABLE)

run: clean all
	clear
	./$(BIN)/$(EXECUTABLE)
	python src/script.py
	@ # cp ./docs/build/main.pdf ./docs.pdf

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	$(CXX) $(CXX_FLAGS) -I $(EIGEN) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*