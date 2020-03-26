all: spglib-1.14.1 eigen-3.3.7

eigen-3.3.7:
	tar -xvf eigen-3.3.7.tar.gz
	touch .

spglib-1.14.1:
	unzip spglib-1.14.1
	mkdir spglib-1.14.1/build
	   (cd ./spglib-1.14.1 && cmake -DCMAKE_INSTALL_PREFIX="")
	   (cd ./spglib-1.14.1 && make)
	   (cd ./spglib-1.14.1 && make DESTDIR=./build install)

clean:
	@$(RM) -r spglib-1.14.1
	@$(RM) -r eigen-3.3.7
