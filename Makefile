CFLAGS=-Wall


all: src/matrix.cc src/align.cc src/editdist.c src/gzstream/libgzstream.a
	@echo "******************************************************************"
	@echo "This is the wrong way to install mabs. The correct way is to say: "
	@echo "# python setup.py install" 
	@echo "******************************************************************"
	@echo 
	python2.2 setup.py build
	cp  build/lib*/*.so .
	#cp build/scripts*/mabs .


modules/matrix.so: src/matrix.cc

modules/editdist.so: src/editdist.c

modules/align.so: src/align.cc src/gzstream/libgzstream.a

src/gzstream/libgzstream.a:
	make -C src/gzstream/

clean:
	rm -rf modules/matrix.so modules/align.so modules/*.pyc *.pyc build
