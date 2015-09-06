# by R.Lie, Nov 01, 2002

include /usr/local/AFEPack/Make.global_options

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
LDFLAGS += -L/usr/local/AFEPack/library/lib -lAFEPack

all : main

%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) -DDEBUG

main : $(object)
	$(CXX) -o $@ $(object) $(LDFLAGS) $(LIBS)
#	$(CXX) -o $@ $(object) $< $(CXXFLAGS)

clean :
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.dx
	-rm -f *~
	-rm -f *.[nes]

.PHONY : default clean
