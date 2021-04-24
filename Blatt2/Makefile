CXX       := g++
CXX_FLAGS := -std=c++17 -ggdb -O3 -ffast-math
PY 		  := python

BIN     := bin
SRC     := src
EIGEN := ~/Downloads/eigen/
INCLUDE   := include

LIBRARIES :=
EXECUTABLE  := main

TEX := tex
BUILD := build
TEX_FILE := main.tex
TEX_OUTPUT := protokoll.pdf
OUTPUT := output

all: run

finish: $(TEX_OUTPUT)
	cp $(TEX)/$(BUILD)/main.pdf $(TEX_OUTPUT)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	clear
	$(CXX) $(CXX_FLAGS) -I $(INCLUDE) -I $(EIGEN) $^ -o $@ $(LIBRARIES)

run: clean $(BIN)/$(EXECUTABLE)
	clear
	./$(BIN)/$(EXECUTABLE)
	
py: run
	$(PY) $(SRC)/*.py

$(TEX_OUTPUT): py
	latexmk \
	-synctex=1 \
	-interaction=nonstopmode \
	-file-line-error \
	-lualatex \
	-outdir=$(TEX)/$(BUILD) \
	$(TEX)/$(TEX_FILE)


clean:
	-rm $(BIN)/*
	-rm $(BUILD)/*
	-rm $(OUTPUT)/*
	-rm $(TEX_OUTPUT)
	