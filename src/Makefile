#
# $Log$
#


CFLAGS=-Wall -g -I/usr/include/python2.2/ -shared -lstdc++


#all: multiAlign.so matrix.so editdist.so align.so gzstream/libgzstream.a


multiAlign.so: multiAlign.cc multiAlign.h
	$(CC) $(CFLAGS)  multiAlign.cc -o multiAlign.so

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



.PHONY:all clean
