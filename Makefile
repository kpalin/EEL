CFLAGS=-Wall


all: src/matrix.cc src/align.cc src/editdist.c src/gzstream/libgzstream.a
	@echo "******************************************************************"
	@echo "This is the wrong way to install EEL. The correct way is to say: "
	@echo "# python setup.py install" 
	@echo "******************************************************************"
	@echo 
	#make -C src/
	python setup.py build
	cp -u  build/lib*/*.so .
	#cp build/scripts*/eel .


debug:src/matrix.cc src/align.cc src/editdist.c src/gzstream/libgzstream.a
	python setup.py debug build
	cp  build/lib*/*.so .

modules/matrix.so: src/matrix.cc

modules/editdist.so: src/editdist.c

modules/align.so: src/align.cc src/gzstream/libgzstream.a

src/gzstream/libgzstream.a:
	make -C src/gzstream/

clean:
	rm -rf matrix.so align.so eellib/*.pyc *.pyc build
