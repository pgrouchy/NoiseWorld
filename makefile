all: EMM_SymComm EMM_island

EMM_SymComm: master.o mathNet.o params.h

	g++ -O3 -Wunused-variable -o EMM_SymComm master.o mathNet.o

EMM_island: island.o mathNet.o params.h

	g++ -O3 -Wunused-variable -o EMM_island island.o mathNet.o

master.o: master.cpp params.h

	g++ -O3 -Wunused-variable -c -o master.o master.cpp

mathNet.o: mathNet.cpp params.h

	g++ -O3 -Wunused-variable -c -o mathNet.o mathNet.cpp

island.o: island.cpp params.h

	g++ -O3 -Wunused-variable -c -o island.o island.cpp

clean:

	rm mathNet.o master.o island.o EMM_SymComm EMM_island
