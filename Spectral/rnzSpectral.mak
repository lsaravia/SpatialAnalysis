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

O = rnzSpectral.o RWFile.o randomizations.o 

L = -lm -ltiff

rnzSpectral: $(O)
	g++ -o rnzSpectral $(O) $(L)

clean:
	rm rnzSpectral $(O)


# DEPENDENCIES

Spec2D.o: Spec2D.cpp

SpectralAux.o: SpectralAux.cpp

randomizations.o: randomizations.cpp Randomizations.h

RWFile.o: RWFile.cpp RWFile.h

rnzSpectral.o: rnzSpectral.cpp 


