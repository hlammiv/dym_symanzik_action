LDLIBS = -lm
#CFLAGS = -O0 -g -pg -std=gnu99
#CPPFLAGS = -O0 -g -pg 
#CXX=g++-8
CXXFLAGS = -O3 -g -pg -fopenmp -std=c++11
CFLAGS = -O3 -std=gnu99 -Wall -fopenmp 
#CFLAGS = -O3 -std=gnu99 -Wall -fsanitize=address
all: dym expectations operators
dym expectations operators: lattice.o group.o
clean:
	${RM} dym-mod-metro dym-mod-metro-savecfg dym expectations operators lattice.o group.o timer.o

dym-mod: dym-mod.cpp lattice.o group.o

dym-mod-metro: dym-mod-metro.cpp lattice.o group.o timer.o

dym-mod-metro-savecfg: dym-mod-metro-savecfg.cpp lattice.o group.o timer.o

dym-mod-metroOG: dym-mod-metroOG.cpp lattice.o group.o timer.o