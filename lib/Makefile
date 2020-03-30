all: spglib-1.14.1 eigen-3.3.7 pugixml-1.10

eigen-3.3.7:
	tar -xvf eigen-3.3.7.tar.gz
	touch .

spglib-1.14.1:
	unzip spglib-1.14.1
	mkdir spglib-1.14.1/build
	   (cd ./spglib-1.14.1 && cmake -DCMAKE_INSTALL_PREFIX="")
	   (cd ./spglib-1.14.1 && make)
	   (cd ./spglib-1.14.1 && make DESTDIR=./build install)

pugixml-1.10:
	tar -xvf pugixml-1.10.tar.gz
	touch .

clean:
	@$(RM) -r spglib-1.14.1
	@$(RM) -r eigen-3.3.7
	@$(RM) -r pugixml-1.10
