SOURCEDIR=../source
OBJFILES=${SOURCEDIR}/edge.o ${SOURCEDIR}/common.o ${SOURCEDIR}/polygon.o ${SOURCEDIR}/triangle.o ${SOURCEDIR}/triangulation.o
OFILES=edge.o common.o polygon.o triangle.o triangulation.o
CXXFLAGS=-g

dlny:	${OBJFILES}
	mv ${SOURCEDIR}/*.o .
	${CXX} ${OFILES} ${SOURCEDIR}/poly_driver.cpp -o dlny.exe
clean:
	rm -rf *.o
