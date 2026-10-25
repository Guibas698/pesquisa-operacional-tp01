CC = g++
CFLAGS1 = -c -m64 -O2 -fPIC -fexceptions -DNDEBUG -DIL_STD -DILOSTRICTPOD
CFLAGS2 = -m64 -O2 -fPIC -fexceptions -DNDEBUG -DIL_STD -DILOSTRICTPOD -Wall
INCLUDES = -I/home/phmacedoaraujo/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include/ -I/home/phmacedoaraujo/opt/ibm/ILOG/CPLEX_Studio1261/concert/include/
OBJS = readFile.o PL.o
SRCS = readFile.c PL.cpp
LIBS = -L/home/phmacedoaraujo/opt/ibm/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_linux/static_pic -lilocplex -lcplex -L/home/phmacedoaraujo/opt/ibm/ILOG/CPLEX_Studio1261/concert/lib/x86-64_linux/static_pic -lconcert -lm -lpthread
EXEC = PL

compila: PL.cpp readFile.c
	${CC} ${CFLAGS1} ${INCLUDES} ${SRCS}
	${CC} ${CFLAGS2} ${INCLUDES} ${OBJS} -o ${EXEC} ${LIBS}

executa: compila
	./PL "../../Implementacao/instancias/Completos e Assimetricos nos Arcos/br17.atsp"

clean:
	rm -f *.o ${EXEC}
