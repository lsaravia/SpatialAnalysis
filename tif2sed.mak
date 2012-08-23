.PHONY : clean

vpath_src=.. ../../fortify ../../randlib/src ../CaNew
vpath %.c    $(vpath_src)
vpath %.cpp  $(vpath_src)
vpath %.hpp  $(vpath_src)
vpath %.h    $(vpath_src)

# The X11 base dir on your system
X11BASE=/usr/X11R6
# Add directories with X11 include files here
X11INCS=-I$(X11BASE)/include
# put X11 required libraries and directories here
X11LIBS=-L$(X11BASE)/lib -lX11

SDLDEFS = -D__XWIN__

I_DIRS=-I../../fortify -I.. -I../../randlib/src -I../CaNew
#P_DEFS=-DGRAPHICS -DPERIODIC_BOUNDARY

#CFLAGS = -O3 -Wall -Ic:/cpp/fortify -Ic:/cpp/canew -DGRAPHICS -DFORTIFY -fexternal-templates 
CXXFLAGS = -g -Wall $(I_DIRS) $(X11INCS)  $(SDLDEFS) $(P_DEFS)

O = RWFile.o tif2sed.o

L = -lm -ltiff

tif2sed : $(O)
	g++ -o tif2sed $(O) $(L)

clean:
	rm tif2sed $(O)


# DEPENDENCIES

RWFile.o: RWFile.cpp RWFile.h

tif2sed.o: tif2sed.cpp
