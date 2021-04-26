CXX       := g++
CXX_FLAGS := -std=c++17 -ggdb -O3 -ffast-math
PY 		  := python
PDF_VIEWER := okular

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

all: $(BIN)/$(EXECUTABLE)

run: $(BIN)/$(EXECUTABLE)
	clear 
	./$(BIN)/$(EXECUTABLE)

show: $(OUTPUT)/*.pdf
	$(PDF_VIEWER) $(OUTPUT)/*.pdf

finish: $(TEX)/$(BUILD)/*.pdf 
	 cp $(TEX)/$(BUILD)/*.pdf $(TEX_OUTPUT)
	 cp $(TEX)/$(BUILD)/*.pdf ../Protokoll_Blatt2_Overhoff_Herdieckerhoff.pdf
	 

$(OUTPUT)/*.txt: $(BIN)/$(EXECUTABLE)
	clear
	./$(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp | $(BIN) $(OUTPUT)
	clear
	$(CXX) $(CXX_FLAGS) -I $(INCLUDE) -I $(EIGEN) $^ -o $@ $(LIBRARIES)

$(OUTPUT)/*.pdf: $(OUTPUT)/*.txt
	$(PY) $(SRC)/*.py

$(TEX)/$(BUILD)/*.pdf: $(TEX)/*.tex $(TEX)/content/*.tex  | $(OUTPUT)/*.pdf $(TEX)/$(BUILD)
	latexmk \
	-synctex=1 \
	-interaction=nonstopmode \
	-file-line-error \
	-lualatex \
	-outdir=$(TEX)/$(BUILD) \
	$(TEX)/$(TEX_FILE)

$(BIN) $(TEX)/$(BUILD) $(OUTPUT): 
	mkdir $@

clean:
	-rm -r $(BIN)
	-rm -r $(TEX)/$(BUILD)
	-rm -r $(OUTPUT)
	