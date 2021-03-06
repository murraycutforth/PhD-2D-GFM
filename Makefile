#
# Makefile for the 2D-GFM project
#
# External dependencies
#	- The header-only Eigen library. Download from eigen.tuxfamily.org
#	- The blitz++ array library. Working on removing this dependency..
#	- My own PhD-Common headers. Download from https://github.com/murraycutforth/PhD-Common
#	- My PhD-Interface-tracking-methods headers and object files
#
# Based on the file at http://hiltmon.com/
#

# Target
TARGET := 2D-GFM

# External header files
EIGEN := -I ~/Libraries/eigen-v3.3.4
COMMON := -I ./../PhD-Common
RS := -I ./../exact_riemann_solver_stiffenedgas
ITM := $(patsubst ./../%, -I ./../%, $(shell find ./../PhD-Interface-tracking-methods/headers/**/* -name '*.hpp' -exec dirname {} \; | sort | uniq)) 

# External .o files
ITMOBJECTS :=  $(shell find ./../PhD-Interface-tracking-methods/objectfiles/**/* -name '*.o')
EXTOBJS = ../exact_riemann_solver_stiffenedgas/exact_RS_stiffenedgas.o $(ITMOBJECTS)

# External link paths
EXT_LINKPATH := -L /home/raid/mcc74/Libraries/NLopt-v2.4.2/lib -lnlopt -lm

# Folders
BUILDDIR := objectfiles
SRCDIR := sourcecode
INCDIRS := headers
INCLIST := -I./headers

# Files
SOURCES := $(shell find $(SRCDIR) -type f -name *.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.o))

# Folder Lists
BUILDLIST := $(BUILDDIR)
EXTINCLIST := $(ITM) $(RS) $(EIGEN) $(HCLFRAMEWORK) $(COMMON) $(INCLIST)

# Shared Compiler Flags
OPLEVEL := -O3
CFLAGS := -Wall -c -fopenmp -std=c++11 $(OPLEVEL) -g
LINKFLAGS := -fopenmp $(OPLEVEL) $(EXT_LINKPATH)

# Linking Step
$(TARGET): $(OBJECTS)
	@mkdir -p $(BUILDLIST)
	@echo "Linking..."
	@echo "Linking $(TARGET) using options: $(LINKFLAGS)"; g++ $^ $(EXTOBJS) $(LINKFLAGS) -o $(TARGET)
	@echo "Success!"

# Compilation Step
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDLIST)
	@echo "Compiling $< using options: $(CFLAGS)"; g++ $(CFLAGS) $(EXTINCLIST) -o $@ $<

clean:
	@echo "Cleaning $(TARGET)..."; rm $(BUILDDIR)/* $(TARGET)
