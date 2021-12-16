include ../../example.mk

CC=mpic++

LDIR =
OPT=

OBJ = main.o

sph_dlb:
sph_dlb_test: OPT += -DTEST_RUN
sph_dlb_test: sph_dlb

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

sph_dlb: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: sph_dlb

run: sph_dlb_test
	mpirun -np 2 ./sph_dlb

.PHONY: clean all run

clean:
	rm -f *.o *~ core sph_dlb

