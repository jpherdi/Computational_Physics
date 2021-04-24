CXX       := g++
CXX_FLAGS := -std=c++17 -ggdb -O3 -ffast-math

BIN     := bin
SRC     := src
EIGEN := ~/Downloads/eigen/
OUTDIR := ./docs/build
TEX_FILE := ./docs/main.tex

LIBRARIES   := #lib
EXECUTABLE  := main

all: $(BIN)/$(EXECUTABLE)

run: clean all
	clear
	./$(BIN)/$(EXECUTABLE)
	python src/script.py

tex: clean clean_tex all
	./$(BIN)/$(EXECUTABLE)
	python src/script.py
	latexmk -synctex=1 -interaction=nonstopmode -file-line-error -lualatex -outdir=$(OUTDIR) $(TEX_FILE)
	clear
	cp $(OUTDIR)/main.pdf ./docs.pdf

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	$(CXX) $(CXX_FLAGS) -I ./lib/ -I $(EIGEN) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*
	
clean_tex:
	-rm $(BIN)/*
	-rm $(OUTDIR)/*