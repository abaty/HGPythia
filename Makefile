ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11

MKDIR_BIN=mkdir -p $(PWD)/bin

SETPYT=export PYTHIA8DATA=../share/Pythia8/xmldoc

all: mkdirBin setpyt bin/HGPYTHIA.exe

mkdirBin:
	$(MKDIR_BIN)
setpyt:
	$(SETPYT)
bin/HGPYTHIA.exe: src/HGPYTHIA.cc ../lib/libpythia8.a
	$(CXX) src/HGPYTHIA.cc ../lib/libpythia8.a -o bin/HGPYTHIA.exe  -I../include -pedantic -W $(CXXFLAGS) -Wshadow -fPIC -L../lib -Wl,-rpath,../lib -lpythia8  -ldl $(ROOT) -I $(PWD)

clean:
	rm -f *~
	rm -f \#*.*#
	rm -f $(PWD)/include/#*.*#
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/#*.*#
	rm -f $(PWD)/src/*~
	rm -f $(PWD)/bash/#*.*#
	rm -f $(PWD)/bash/*~
	rm -f $(PWD)/paths/#*.*#
	rm -f $(PWD)/paths/*~
	rm -f $(PWD)/bin/*.exe
	rmdir bin
.PHONY: all
