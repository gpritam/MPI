MPICXX := mpic++
IDIR   := include
LDIR   := lib
ODIR   := build
SRCDIR := src
DOCDIR := documents

TARGET := DendriticGrowth.cpp
Nproc  := 4

CXXFLAGS := -I${IDIR}

_DEPENDENCIES := General.h General_MPI_2D.h
DEPENDENCIES  := ${patsubst %,${IDIR}/%,${_DEPENDENCIES}}

_OBJECTFILES := ${patsubst %.h,%.o,${_DEPENDENCIES}} 
_OBJECTFILES += ${patsubst %.cpp,%.o,${TARGET}}
OBJECTFILES  := ${patsubst %,${ODIR}/%,${_OBJECTFILES}}

all: ${OBJECTFILES}
	${MPICXX} -std=c++17 -O3 -o run $^

${ODIR}/%.o: %.cpp ${DEPENDENCIES}
	${MPICXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

${ODIR}/%.o: ${SRCDIR}/%.cpp ${DEPENDENCIES}
	${MPICXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

run:
	mpirun -np ${Nproc} ./run 2>/dev/null

.PHONY: clean run document

document:
	@xelatex -synctex=1 -interaction=nonstopmode ${DOCDIR}/Document.tex >> /dev/null
	@rm -rf ${DOCDIR}/*.aux ${DOCDIR}/*.log ${DOCDIR}/*.gz ${DOCDIR}/*.out
	@rm -rf *.aux *.log *.gz *.out
	@mv *.pdf ${DOCDIR}/

clean:
	rm -rf ${ODIR}/*.o ${SRCDIR}/*~ *~ run ${IDIR}/*~ Output/*.tec Output/*.vtk Output/*.visit
