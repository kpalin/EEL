#
# $Log$
# Revision 1.2  2004/07/30 12:17:49  kpalin
# Added alnColumn and usage of top level makefile.
#
# Revision 1.1  2004/07/12 11:23:37  kpalin
# Multiple alignment added.
#
#


CFLAGS=-Wall -g -I/usr/include/python2.2/ -shared -lstdc++
BASEDIR=$(PWD)/../

debug:
	cd $(BASEDIR); \
	python setup.py  debug install --prefix=/home/kpalin/tyot/comparative/EEL/subzero/

all:
	cd $(BASEDIR);\
	python setup.py debug install --prefix=/home/kpalin/fs/EEL/

#all: multiAlign.so matrix.so editdist.so align.so gzstream/libgzstream.a


multiAlign.so: multiAlign.cc multiAlign.h alignedCols.cc alignedCols.h
	gcc -pthread -fno-strict-aliasing -DNDEBUG -O2 -g -pipe -march=i386 -mcpu=i686 -D_GNU_SOURCE -fPIC -fPIC -I/usr/include/python2.3 -c multiAlign.cc -DSAVE_MEM -DSAVE_MEM_LIMIT=536870912 -Isrc/gzstream -DHAVE_GZSTREAM=1 -Wall -O3 -Wall -shared
#	make -C `pwd`/../ debug
#	$(CC) $(CFLAGS)  multiAlign.cc -o multiAlign.so

matrix.so: matrix.cc
	$(CC) $(CFLAGS)  $? -o $@

editdist.so: editdist.c
	$(CC) $(CFLAGS)  $? -o $@

align.so: align.cc gzstream/libgzstream.a
	$(CC) $(CFLAGS)  $? -o $@

gzstream/libgzstream.a:
	make -C src/gzstream/

clean:
	rm -f matrix.so editdist.so align.so gzstream/libgzstream.a
	rm -r $(BASEDIR)/build



.PHONY:all debug clean
