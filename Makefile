################################################################################
#
# Makefile for the pseudo-spectral code TURBID
#
################################################################################



###
# Set the lib path for the libraries of FFTW, Armadillo and P3DFFT
###
WORKPATH    := .
BUILDPATH   := $(WORKPATH)/build
BINPATH     := $(WORKPATH)/bin
SOURCEPATH  := $(WORKPATH)/src
LIBPATH     := $(WORKPATH)/lib
FFTWPATH    := Y_PATH/Library/fftw/3.3.6-pl2
ARMAPATH    := Y_PATH/Library/armadillo/7.950.1
P3DFFTPATH  := Y_PATH/Library/p3dfft/2.7.5



###
#
###
CC          := mpic++
CFLAGS      := -std=c++11 -O3 -fopenmp -Wall -DHAVE_CONFIG_H \
               -DGNU -DOPENMP -DONED -DMEASURE -DSTRIDE1 -DFFTW \
               -g -DARMA_DONT_USE_WRAPPER -DUSE_MKL
LIBS        := -larmadillo -lpseudospectral \
              -lm -fopenmp -lp3dfft -lfftw3 -lgfortran -lmpi_mpifh
CPPFLAGS    := -I$(WORKPATH)/include \
               -I$(FFTWPATH)/include \
               -I$(ARMAPATH)/usr/include \
               -I$(P3DFFTPATH)/include \
               -I$(SOURCEPATH)
LDFLAGS     := -L$(FFTWPATH)/lib \
	       -L$(ARMAPATH)/usr/lib64 \
               -L$(P3DFFTPATH)/lib \
               -L$(WORKPATH)/lib



###
#
###
SOURCEEXT   := cc
SOURCEFILE  := $(shell find $(SOURCEPATH)/ -type f -name '*.$(SOURCEEXT)')
SOURCEHEAD  := $(shell find $(SOURCEPATH)/ -type f -name '*.h')
OBJECTS     := $(patsubst $(SOURCEPATH)/%,$(BUILDPATH)/%,$(SOURCEFILE:.$(SOURCEEXT)=.o))
TARGETLIB   := $(WORKPATH)/lib/libpseudospectral.so




###
#
###
.DEFAULT: test



###
#
###
.PHONY  : lib
lib     : $(TARGETLIB)

$(TARGETLIB) : $(OBJECTS)
	$(CC) -shared -o $@ $^

$(OBJECTS) : $(BUILDPATH)/%.o : $(SOURCEPATH)/%.$(SOURCEEXT)
	$(CC) $(CFLAGS) -fPIC -c $< $(CPPFLAGS) -o $@



###
#
###
.PHONY: flat
flat : $(BINPATH)/channelflat

$(BINPATH)/channelflat : $(BUILDPATH)/channelflat.o lib
	$(CC) $(BUILDPATH)/channelflat.o -o $(BINPATH)/channelflat $(LDFLAGS) $(LIBS)

$(BUILDPATH)/channelflat.o : $(WORKPATH)/main/channelflat.$(SOURCEEXT)
	$(CC) -c $(WORKPATH)/main/channelflat.$(SOURCEEXT) $(CFLAGS) -o $(BUILDPATH)/channelflat.o $(CPPFLAGS)



###
#
###
.PHONY: clean
clean:
	rm -f $(BUILDPATH)/* &
	rm -f $(LIBPATH)/* &
	rm -f $(BINPATH)/*




