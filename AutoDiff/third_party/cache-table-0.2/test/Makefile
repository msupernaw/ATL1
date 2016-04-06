

bins = map_unittest 
sources = map_unittest.cpp

#################################################

FLAGS = -O3 -Wall -fomit-frame-pointer -DNDEBUG 
INCLUDES = -I..
LIBS = -lstdc++

CXXFLAGS= $(FLAGS) $(INCLUDES) 
CXX=g++

# sources = $(bins:=.cpp)

all: $(bins)

$(bins): $(sources)
	$(CXX) $(FLAGS) $(INCLUDES) $(LIBS)  $< -o $@ 

clean: 
	rm -rf *~ *.d *.o $(bins)

#%.d: %.cpp
#	set -e; $(CXX) $(INCLUDES) -MM $(CPPFLAGS) $< \
#		| sed 's/\($*\)\.o[ :]*/\1 $@ : /g' > $@; \
#	        [ -s $@ ] || rm -f $@
#
# include $(sources:.cpp=.d)
