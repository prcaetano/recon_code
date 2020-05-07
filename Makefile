CC     = g++
OPT    = -funroll-loops -O -fopenmp -fPIC -DSKIPRAW # -DREADWEIGHT


all: recon.so recon


recon: recon.o io.o multigrid.o grid.o shift.o smooth.o
	$(CC) $(OPT) recon.o io.o multigrid.o grid.o shift.o smooth.o -o recon -lfftw3_omp -lfftw3 -lm
recon.so: recon.o io.o multigrid.o grid.o shift.o smooth.o
	$(CC) -shared $^ -o $@ -Wl,--whole-archive -lfftw3_omp -lfftw3 -lm -Wl,--no-whole-archive
tex: notes.tex notes.bib
	pdflatex notes
	bibtex   notes
	pdflatex notes
	pdflatex notes
	rm -f    notes.toc



recon.o: recon.cpp global.h lcdm.h multigrid.h
	$(CC) $(OPT) -c recon.cpp
io.o: io.cpp global.h lcdm.h
	$(CC) $(OPT) -c io.cpp
multigrid.o: multigrid.cpp global.h lcdm.h multigrid.h
	$(CC) $(OPT) -c multigrid.cpp
grid.o: grid.cpp global.h lcdm.h
	$(CC) $(OPT) -c grid.cpp
shift.o: shift.cpp global.h lcdm.h
	$(CC) $(OPT) -c shift.cpp
smooth.o: smooth.cpp global.h
	$(CC) $(OPT) -I$(FFTW_INC) -L$(FFTW_DIR) -c smooth.cpp


.(PHONY) clean:
	rm -f notes.aux notes.log notes.out notes.toc *.o
