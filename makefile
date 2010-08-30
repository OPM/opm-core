PEDANTIC=-W -Wall -pedantic # -W -Wformat-nonliteral \
        #-Wcast-align -Wpointer-arith -Wbad-function-cast \
        #-Wmissing-prototypes -Wstrict-prototypes \
        #-Wmissing-declarations -Winline -Wundef -Wnested-externs\
        #-Wcast-qual -Wshadow -Wconversion -Wwrite-strings\
        #-Wno-conversion -Wchar-subscripts -Wredundant-decls\



GCC     = gcc -fPIC 
CC      = $(GCC)
INC     = -I.

#rememeber export LD_LIBRARY_PATH=/home/jrn/devel/c/repo/preprocess
LIB    =  -L.
CFLAGS = -g $(INC) -Wall $(PEDANTIC)
SHELL  = /bin/bash
DEFAULT_LIB_INSTALL_PATH = /home/jrn/devel/c/repo/opm-gridprocessing
all:    tags depend libpreprocess                                        

OBJ   = preprocess.o uniquepoints.o facetopology.o sparsetable.o newinterface.o ../reorder-utils/grid.o

libpreprocess:	libpreprocess.so.1.0.1
libpreprocess.so.1.0.1: $(OBJ)
	$(CC) -shared -Wl,-soname,libpreprocess.so.1,\
	-rpath,$(DEFAULT_LIB_INSTALL_PATH) -o  $@  $(OBJ)  -lc $(LIB)
	ln -s libpreprocess.so.1.0.1 libpreprocess.so.1
	ln -s libpreprocess.so.1     libpreprocess.so
.PHONY: clean depend all

clean:; rm -f *~ $(OBJ) libpreprocess.so.1.0.1 libpreprocess.so.1 \
          libpreprocess.so TAGS; makedepend

tags : $(OBJ:.o=.c)  $(wildcard *.h)
	etags *.c *.h

depend :; makedepend $(INC) -f makefile *.c 2>/dev/null

# DO NOT DELETE THIS LINE -- make depend depends on it.                         
