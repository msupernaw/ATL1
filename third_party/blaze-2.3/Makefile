#==================================================================================================
#
#  Makefile for the Blaze library
#
#  Copyright (C) 2013 Klaus Iglberger - All Rights Reserved
#
#  This file is part of the Blaze library. You can redistribute it and/or modify it under
#  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
#  forms, with or without modification, are permitted provided that the following conditions
#  are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this list of
#     conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice, this list
#     of conditions and the following disclaimer in the documentation and/or other materials
#     provided with the distribution.
#  3. Neither the names of the Blaze development group nor the names of its contributors
#     may be used to endorse or promote products derived from this software without specific
#     prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
#  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
#  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
#  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
#  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
#  DAMAGE.
#
#==================================================================================================


# Compiler configuration
export CXX      = g++
export CXXFLAGS = -Werror -Wall -Wextra -Wshadow -Woverloaded-virtual -ansi -O3 -DNDEBUG

# Exporting the installation path and object path
export INSTALL_PATH = /home/matthew/NetBeansProjects/ATL/ATL/third_party/blaze-2.3
export OBJECT_PATH  = /home/matthew/NetBeansProjects/ATL/ATL/third_party/blaze-2.3/obj

# Exporting the boost include path
export BOOST_INCLUDE_PATH = 


# Build rules
default: static
	@echo

static: util math
	@echo
	@echo "Creating static library"
	@ar crus $(INSTALL_PATH)/lib/libblaze.a $(OBJECT_PATH)/*.o
	@echo "... finished!"

shared: util math
	@echo
	@echo "Creating shared library"
	@g++ -shared -o $(INSTALL_PATH)/lib/libblaze.so $(OBJECT_PATH)/*.o  -lboost_system -lboost_thread
	@echo "... finished!"

util:
	@echo
	@echo "Setup of the utility module..."
	@$(MAKE) --no-print-directory -C $(INSTALL_PATH)/src/util
	@echo "...Setup of the utility module complete!"

math:
	@echo
	@echo "Setup of the math module..."
	@$(MAKE) --no-print-directory -C $(INSTALL_PATH)/src/math
	@echo "...Setup of the math module complete!"


# Clean up rules
clean:
	@echo "Cleaning up..."
	@$(MAKE) --no-print-directory -C $(INSTALL_PATH)/src/util clean
	@$(MAKE) --no-print-directory -C $(INSTALL_PATH)/src/math clean
	@$(RM) $(INSTALL_PATH)/lib/libblaze.a $(INSTALL_PATH)/lib/libblaze.so $(OBJECT_PATH)/*.o


# Doxygen documentation
doc:
	@doxygen Doxyfile

doc_clean:
	@$(RM) -r $(INSTALL_PATH)/doc/html/ $(INSTALL_PATH)/doc/latex/


# Setting the independent commands
.PHONY: default clean doc doc_clean
