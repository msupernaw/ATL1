## -*- Makefile -*-
##
## User: matthew
## Time: Jun 25, 2015 9:07:06 AM
## Makefile created by Oracle Solaris Studio.
##
## This file is generated automatically.
##


#### Compiler and tool definitions shared by all build targets #####
CCC = g++
CXX = g++
BASICOPTS = -O3 -s -m64
CCFLAGS = -std=c++11
CXXFLAGS =-O3 -std=c++11
CCADMIN = 


# Define the target directories.
TARGETDIR_autodifftests=build


all: $(TARGETDIR_autodifftests)/autodifftests

## Target: autodifftests
OBJS_autodifftests =  \
	$(TARGETDIR_autodifftests)/main.o
USERLIBS_autodifftests = $(SYSLIBS_autodifftests) 
DEPLIBS_autodifftests =  
LDLIBS_autodifftests = $(USERLIBS_autodifftests)


# Link or archive
$(TARGETDIR_autodifftests)/autodifftests: $(TARGETDIR_autodifftests) $(OBJS_autodifftests) $(DEPLIBS_autodifftests)
	$(LINK.cc) $(CCFLAGS_autodifftests) $(CPPFLAGS_autodifftests) -o $@ $(OBJS_autodifftests) $(LDLIBS_autodifftests)


# Compile source files into .o files
$(TARGETDIR_autodifftests)/main.o: $(TARGETDIR_autodifftests) main.cpp
	$(COMPILE.cc) $(CCFLAGS_autodifftests) $(CPPFLAGS_autodifftests) -o $@ main.cpp



#### Clean target deletes all generated files ####
clean:
	rm -f \
		$(TARGETDIR_autodifftests)/autodifftests \
		$(TARGETDIR_autodifftests)/main.o
	$(CCADMIN)
	rm -f -r $(TARGETDIR_autodifftests)


# Create the target directory (if needed)
$(TARGETDIR_autodifftests):
	mkdir -p $(TARGETDIR_autodifftests)


# Enable dependency checking
.KEEP_STATE:
.KEEP_STATE_FILE:.make.state.GNU-amd64-Linux

run : ${Target}
	./$(TARGETDIR_autodifftests)/autodifftests

