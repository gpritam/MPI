MPICXX := mpic++
IDIR   := include
LDIR   := lib
ODIR   := build
SRCDIR := src
DOCDIR := documents

CXXFLAGS := -I${IDIR}
INCPATH = -I/home/pritam/Hypre/hypre-2.11.2/hypre/include/
LIBPATH = -L/home/pritam/Hypre/hypre-2.11.2/hypre/lib/ -lHYPRE -lm

TARGET := TestHelmholtz1D_MPI.cpp
Nproc  := 4

_DEPENDENCIES := General.h HypreSolvers.h 
_DEPENDENCIES += ${patsubst Test%.cpp,%.h,${TARGET}}
DEPENDENCIES  := ${patsubst %,${IDIR}/%,${_DEPENDENCIES}}

_OBJECTFILES := ${patsubst %.h,%.o,${_DEPENDENCIES}} ${patsubst %.cpp,%.o,${TARGET}}
OBJECTFILES  := ${patsubst %,${ODIR}/%,${_OBJECTFILES}}

all: ${OBJECTFILES}
	${MPICXX} -std=c++17 -O3 -o run $^ ${LIBPATH}

${ODIR}/%.o: ${SRCDIR}/%.cpp ${DEPENDENCIES}
	${MPICXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

${ODIR}/${patsubst Test%.cpp,%.o,${TARGET}}: ${SRCDIR}/${patsubst Test%,%,${TARGET}} ${DEPENDENCIES}
	${MPICXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS} ${INCPATH} ${LIBPATH}

${ODIR}/HypreSolvers.o: ${SRCDIR}/HypreSolvers.cpp ${IDIR}/HypreSolvers.h ${IDIR}/General.h
	${MPICXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS} ${INCPATH} ${LIBPATH}

${ODIR}/${patsubst %.cpp,%.o,${TARGET}}: ${TARGET} ${DEPENDENCIES}
	${MPICXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS} ${INCPATH} ${LIBPATH}

run:
	mpirun -np ${Nproc} ./run 2>/dev/null

.PHONY: clean document run

document:
	@xelatex -synctex=1 -interaction=nonstopmode ${DOCDIR}/Document.tex >> /dev/null
	@rm -rf ${DOCDIR}/*.aux ${DOCDIR}/*.log ${DOCDIR}/*.gz ${DOCDIR}/*.out
	@rm -rf *.aux *.log *.gz *.out
	@mv *.pdf ${DOCDIR}/

clean:
	rm -rf ${ODIR}/*.o ${SRCDIR}/*~ *~ run ${IDIR}/*~ Output/*.tec Output/*.vtk Output/*.visit
