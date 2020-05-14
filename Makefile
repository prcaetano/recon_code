CC     = g++
OPT    = -funroll-loops -O3 -fopenmp -fPIC -ffast-math -DSKIPRAW -I./include #-DDISTANT_OBSERVER_ZAXIS # -DREADWEIGHT

all: ./lib/recon.so recon


recon: ./lib/recon.o ./lib/io.o ./lib/multigrid.o ./lib/grid.o ./lib/shift.o ./lib/smooth.o
	$(CC) $(OPT) ./lib/recon.o ./lib/io.o ./lib/multigrid.o ./lib/grid.o ./lib/shift.o ./lib/smooth.o -o recon -lfftw3_omp -lfftw3 -lm
./lib/recon.so: ./lib/recon.o ./lib/io.o ./lib/multigrid.o ./lib/grid.o ./lib/shift.o ./lib/smooth.o
	$(CC) -shared $^ -o $@ -Wl,--whole-archive -lfftw3_omp -lfftw3 -lm -Wl,--no-whole-archive
tex: ./tex/notes.tex ./tex/notes.bib
	pdflatex ./tex/notes
	bibtex   ./tex/notes
	pdflatex ./tex/notes
	pdflatex ./tex/notes
	rm -f    ./tex/notes.toc



./lib/recon.o: ./src/recon.cpp ./include/global.h ./include/lcdm.h ./include/multigrid.h
	$(CC) $(OPT) -c ./src/recon.cpp -o ./lib/recon.o
./lib/io.o: ./src/io.cpp ./include/global.h ./include/lcdm.h
	$(CC) $(OPT) -c ./src/io.cpp -o ./lib/io.o
./lib/multigrid.o: ./src/multigrid.cpp ./include/global.h ./include/lcdm.h ./include/multigrid.h
	$(CC) $(OPT) -c ./src/multigrid.cpp -o ./lib/multigrid.o
./lib/grid.o: ./src/grid.cpp ./include/global.h ./include/lcdm.h
	$(CC) $(OPT) -c ./src/grid.cpp -o ./lib/grid.o
./lib/shift.o: ./src/shift.cpp ./include/global.h ./include/lcdm.h
	$(CC) $(OPT) -c ./src/shift.cpp -o ./lib/shift.o
./lib/smooth.o: ./src/smooth.cpp ./include/global.h
	$(CC) $(OPT) -I$(FFTW_INC) -L$(FFTW_DIR) -c ./src/smooth.cpp -o ./lib/smooth.o


.(PHONY) clean:
	rm -f notes.aux notes.log notes.out notes.toc ./lib/*.o
