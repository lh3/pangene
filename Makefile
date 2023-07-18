CC=			gcc
CXX=		g++
CFLAGS=		-std=c99 -g -Wall -O2
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		sys.o dict.o option.o format.o read.o hit.o vertex.o branch.o graph.o
PROG=		pangene
LIBS=		-lpthread -lz -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

pangene:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

branch.o: pgpriv.h pangene.h
dict.o: pgpriv.h pangene.h khashl.h
format.o: pgpriv.h pangene.h
graph.o: pgpriv.h pangene.h ksort.h
hit.o: pgpriv.h pangene.h ksort.h
main.o: pgpriv.h pangene.h ketopt.h
option.o: pangene.h
read.o: pgpriv.h pangene.h kseq.h
sys.o: pgpriv.h pangene.h
vertex.o: pgpriv.h pangene.h
