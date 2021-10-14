CC     = CC
OPT    = -funroll-loops -O -fopenmp -I./include -DSKIPRAW -DWRITECARTESIAN -DPRINT_GRID # -DREADWEIGHT #-DWRITECARTESIAN



recon: recon.o io.o multigrid.o grid.o shift.o smooth.o
	$(CC) $(OPT) recon.o io.o multigrid.o grid.o shift.o smooth.o -o recon -lfftw3_omp -lfftw3 -lm

recon.o: ./src/recon.cpp ./include/global.h ./include/lcdm.h ./include/multigrid.h
	$(CC) $(OPT) -c ./src/recon.cpp
io.o: ./src/io.cpp ./include/global.h ./include/lcdm.h
	$(CC) $(OPT) -c ./src/io.cpp
multigrid.o: ./src/multigrid.cpp ./include/global.h ./include/lcdm.h ./include/multigrid.h
	$(CC) $(OPT) -c ./src/multigrid.cpp
grid.o: ./src/grid.cpp ./include/global.h ./include/lcdm.h
	$(CC) $(OPT) -c ./src/grid.cpp
shift.o: ./src/shift.cpp ./include/global.h ./include/lcdm.h
	$(CC) $(OPT) -c ./src/shift.cpp
smooth.o: ./src/smooth.cpp ./include/global.h
	$(CC) $(OPT) -I$(FFTW_INC) -L$(FFTW_DIR) -c ./src/smooth.cpp



tex:	./tex/notes.tex ./tex/notes.bib
	pdflatex ./tex/notes
	bibtex   ./tex/notes
	pdflatex ./tex/notes
	pdflatex ./tex/notes
	rm -f    ./tex/notes.toc


all:	recon tex


.(PHONY) clean:
	rm -f ./tex/notes.aux ./tex/notes.log ./tex/notes.out ./tex/notes.toc *.o
