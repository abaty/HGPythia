ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11

MKDIR_BIN=mkdir -p $(PWD)/bin

SETPYT=export PYTHIA8DATA=/Users/austinbaty/Desktop/HG_PYTHIA/pythia8243/share/Pythia8/xmldoc

all: mkdirBin setpyt bin/HGPYTHIA.exe

mkdirBin:
	$(MKDIR_BIN)
setpyt:
	$(SETPYT)
bin/HGPYTHIA.exe: src/HGPYTHIA.cc /Users/austinbaty/Desktop/HG_PYTHIA/pythia8243/lib/libpythia8.a
	$(CXX) src/HGPYTHIA.cc /Users/austinbaty/Desktop/HG_PYTHIA/pythia8243/lib/libpythia8.a -o bin/HGPYTHIA.exe  -I/Users/austinbaty/Desktop/HG_PYTHIA/pythia8243/include -pedantic -W $(CXXFLAGS) -Wshadow -fPIC -L/Users/austinbaty/Desktop/HG_PYTHIA/pythia8243/lib -Wl,-rpath,/Users/austinbaty/Desktop/HG_PYTHIA/pythia8243/lib -lpythia8  -ldl $(ROOT) -I $(PWD)

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
