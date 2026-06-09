LDLIBS = -lm
CXXFLAGS = -O3 -g -pg -fopenmp -std=c++11
CFLAGS   = -O3 -std=gnu99 -Wall -fopenmp

all: dym-mod-metro

dym-mod-metro: dym-mod-metro.cpp lattice.o group.o timer.o
dym-mod-metro-savecfg: dym-mod-metro-savecfg.cpp lattice.o group.o timer.o
dym-mod-metroOG: dym-mod-metroOG.cpp lattice.o group.o timer.o
llr: llr.cpp lattice.o group.o timer.o

# Group-file format checker (uses the real load_group from group.o)
verify_group: verify_group.c group.o
	$(CC) $(CFLAGS) $< group.o $(LDLIBS) -o $@

clean:
	${RM} dym-mod-metro dym-mod-metro-savecfg dym-mod-metroOG verify_group \
	      lattice.o group.o timer.o
