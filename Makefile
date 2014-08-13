#!/bin/make
# GNUmakefile using implicit rules and standard definitions
SHELL=/bin/sh

# Uncomment the following line for root functionality
USEROOT = 1
# Uncomment this line if processing Rev. D data
#REVISIOND = 1
# Uncomment this line if processing Rev. F data
REVISIONF = 1
# Uncomment this line for a more verbose scan
# CXXFLAGS += -DVERBOSE
# Undefine to make a "online" version
# ONLINE = 1 

# Define to use Gamma-Gamma gates in GeProcessor
# This turns on Gamma-Gamma angular distribution
# and Gamma-Gamma-Gamma gates
# GGATES = 1

# Define to see debugging information for TreeCorrelator
#DEBUG = 1

# Use gfortran
ifeq ($(HHIRF_GFORTRAN),)
HHIRF_GFORTRAN = 1
endif

#These will set the analysis used on the waveforms
#Uncomment this line to use the Pulse Fitting routine
PULSEFIT = 1
#Uncomment this line to use the cfd
#DCFD = 1


#------- instruct make to search through these
#------- directories to find files
vpath %.f scan/ 
vpath %.hpp include/
vpath %.h include/
vpath %.icc include/
vpath %.cpp src/
vpath %.o obj/

ifeq ($(HHIRF_DIR),)
HHIRF_DIR = /usr/hhirf-intel64
endif

ifneq ($(LIBS_IN_HHIRF),)
ACQ2_LIBDIR = $(HHIRF_DIR)
else
ifeq ($(ACQ2_LIBDIR),)
ifneq ($(ACQ2_DIR),) 
ACQ2_LIBDIR = $(ACQ2_DIR)
else
ACQ2_LIBDIR = /usr/acq2n/lib
endif
endif
endif

LIBS = $(HHIRF_DIR)/scanorlib.a $(HHIRF_DIR)/orphlib.a \
       $(ACQ2_LIBDIR)/acqlib.a  $(ACQ2_LIBDIR)/ipclib.a

OutPutOpt     = -o # keep whitespace after "-o"

#------- define file suffixes
fSrcSuf   = f
cSrcSuf   = c
c++SrcSuf = cpp
cxxSrcSuf = cxx

#------- define compilers
#define to compile with gfortran (>=4.2) if required for the hhirf libs
ifneq ($(HHIRF_GFORTRAN),1)
FC        = g77
else
FC        = gfortran
endif

GCC       = gcc 
CXX       = g++
LINK.o    = $(FC) $(LDFLAGS)

# -Dnewreadout is needed to account for a change to pixie16 readout
# structure change on 03/20/08.  Remove for backwards compatability
#
# for debug and profiling add options -g -pg
# and remove -O
#------- define basic compiler flags, no warnings on code that is not our own
#FFLAGS   += -O3
FFLAGS += -g
GCCFLAGS += -fPIC $(CINCLUDEDIRS) -Dnewreadout
CXXFLAGS += -Wall -g -fPIC $(CINCLUDEDIRS) -Dnewreadout

ifdef REVISIOND
CXXFLAGS += -DREVD
endif
ifdef REVISIONF
CXXFLAGS += -DREVF
endif
ifdef ONLINE
CXXFLAGS += -DONLINE
endif

#------- include directories for the pixie c files
CINCLUDEDIRS  = -Iinclude

#------- basic linking instructions
LDLIBS   += -lm -lstdc++
ifdef PULSEFIT
LDLIBS   += -lgsl -lgslcblas
CXXFLAGS += -Dpulsefit
else ifdef DCFD
CXXFLAGS += -Ddcfd
endif

ifeq ($(FC),gfortran)
FFLAGS	 += -fsecond-underscore
LDLIBS	 += -lgfortran
#GCCFLAGS += -O3
#CXXFLAGS += -O3 -DLINK_GFORTRAN
CXXFLAGS += -DLINK_GFORTRAN -g
else
LDFLAGS += -g77libs
LDLIBS	+= -lg2c
endif

#-------- define file variables -----------------------
# objects from fortran
MESSLOGO         = messlog.o
MILDATIMO        = mildatim.o
SCANORO          = scanor.o
SET2CCO          = set2cc.o

# objects from cpp
#General Objects
CHANEVENTO       = ChanEvent.o
CHANIDENTIFIERO  = ChanIdentifier.o
CORRELATORO      = Correlator.o
DETECTORDRIVERO  = DetectorDriver.o
DETECTORLIBRARYO = DetectorLibrary.o
DETECTORSUMMARYO = DetectorSummary.o
EVENTPROCESSORO  = EventProcessor.o
INITIALIZEO      = Initialize.o
MAPFILEO         = MapFile.o
PATHHOLDERO	 = PathHolder.o
PIXIEO           = PixieStd.o
PLACESO          = Places.o
PLACEBUILDERO    = PlaceBuilder.o
PUGIXMLO	 = pugixml.o
RANDOMPOOLO      = RandomPool.o
RAWEVENTO        = RawEvent.o
STATSDATAO       = StatsData.o
TIMINGINFOO      = TimingInformation.o
TREECORRELATORO  = TreeCorrelator.o

#Plot Objects
PLOTSREGISTERO   = PlotsRegister.o
PLOTSO           = Plots.o

#Trace Related Objects
CFDANALYZERO     = CfdAnalyzer.o
DOUBLETRACEO     = DoubleTraceAnalyzer.o
FITTINGANALYZERO = FittingAnalyzer.o
TAUANALYZERO     = TauAnalyzer.o
TRACEO		 = Trace.o
TRACEEXTRACTORO  = TraceExtractor.o
TRACEFILTERO     = TraceFilterer.o
TRACEPLOTO       = TracePlotter.o
TRACESUBO        = TraceAnalyzer.o
WAVEFORMSUBO     = WaveformAnalyzer.o

#Dectector Processor Objects
BETAPROCESSORO   = BetaProcessor.o
DSSDPROCESSORO   = DssdProcessor.o
GEPROCESSORO     = GeProcessor.o
IONCHAMBERPROCESSORO = IonChamberProcessor.o
ISSDPROCESSORO   = ImplantSsdProcessor.o
LIQUIDPROCESSORO = LiquidProcessor.o
LOGICPROCESSORO  = LogicProcessor.o
MCPPROCESSORO    = McpProcessor.o
MTCPROCESSORO    = MtcProcessor.o
NEUTRONPROCESSORO = NeutronProcessor.o
POSITIONPROCESSORO = PositionProcessor.o
PULSERPROCESSORO = PulserProcessor.o
SSDPROCESSORO    = SsdProcessor.o
VALIDPROCESSORO  = ValidProcessor.o
VANDLEPROCESSORO = VandleProcessor.o

#ROOT Objects
#SCINTROOTO       = ScintROOT.o
ROOTPROCESSORO   = RootProcessor.o
VANDLEROOTO      = VandleROOT.o

ifdef REVISIOND
READBUFFDATAO    = ReadBuffData.RevD.o
else
ifdef REVISIONF
READBUFFDATAO    = ReadBuffData.RevD.o
else
READBUFFDATAO    = ReadBuffData.o
endif
endif

#----- list of objects
#Objects from Fortran
FORT_OBJS = $(MESSLOGO) $(MILDATIMO) $(SCANORO) $(SET2CCO) 
#General Objects
CXX_OBJS += $(CHANEVENTO) $(CHANIDENTIFIERO) $(CORRELATORO) $(DETECTORDRIVERO) \
	$(DETECTORLIBRARYO) $(DETECTORSUMMARYO) $(EVENTPROCESSORO) $(INITIALIZEO) \
	$(MAPFILEO) $(PATHHOLDERO) $(PIXIEO) $(PLACESO) $(PLACEBUILDERO) $(PUGIXMLO) \
	$(RANDOMPOOLO) $(RAWEVENTO) $(READBUFFDATAO) $(STATSDATAO) \
	$(TIMINGINFOO) $(TREECORRELATORO)
#Plot Objects
CXX_OBJS += $(HISTOGRAMMERO) $(PLOTSREGISTERO) $(PLOTSO)
#Trace Objects
CXX_OBJS += $(DOUBLETRACEO) $(TAUANALYZERO) $(TRACEO) $(TRACEEXTRACTORO) \
	$(TRACEFILTERO) $(TRACEPLOTO) $(TRACESUBO) $(WAVEFORMSUBO) 
#Detector Objects
CXX_OBJS += $(BETAPROCESSORO) $(DSSDPROCESSORO) $(GEPROCESSORO) \
	$(IONCHAMBERPROCESSORO) $(ISSDPROCESSORO) $(LIQUIDPROCESSORO) \
	$(LOGICPROCESSORO) $(MCPPROCESSORO) $(MTCPROCESSORO) \
	$(NEUTRONPROCESSORO) $(POSITIONPROCESSORO) $(PULSERPROCESSORO) \
	$(SSDPROCESSORO) $(VALIDPROCESSORO) $(VANDLEPROCESSORO)

ifdef PULSEFIT
CXX_OBJS += $(FITTINGANALYZERO)
else ifdef DCFD
CXX_OBJS += $(CFDANALYZERO) 
endif

#---------- Change the executable name if necessary
ifdef ONLINE
PIXIE = pixie_ldf_c_online$(ExeSuf)
else
PIXIE = pixie_ldf_c$(ExeSuf)
endif

#---------- Adjust compilation if ROOT capability is desired
ifdef USEROOT
CXX_OBJS  += $(ROOTPROCESSORO) $(VANDLEROOTO) $(SCINTROOTO)
PIXIE = pixie_ldf_c_root$(ExeSuf)
ROOTCONFIG   := root-config

#no uncomment ROOTCLFAGS   := $(filter-out pthread,$(ROOTCFLAGS))
CXXFLAGS     += $(shell $(ROOTCONFIG) --cflags) -Duseroot
LDFLAGS      += $(shell $(ROOTCONFIG) --ldflags)
LDLIBS += $(shell $(ROOTCONFIG) --libs)
#LDLIBS       += -L/usr/local/root/lib $(shell $(ROOTCONFIG) --libs)
#LDLIBS		+= -L/usr/local/root/lib -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lMathCore -lThread -pthread -lm -ldl -rdynamic
endif

#------------ Compile with Gamma-Gamma gates support in GeProcessor
ifdef GGATES
CXXFLAGS	+= -DGGATES
endif

#------------ Compile with debug information for TreeCorrelator
ifdef GGATES
CXXFLAGS	+= -DDEBUG
endif

#---------- Update some information about the object files 
FORT_OBJDIR = obj/fortran
FORT_OBJS_W_DIR = $(addprefix $(FORT_OBJDIR)/,$(FORT_OBJS))
CXX_OBJDIR = obj/c++
CXX_OBJS_W_DIR = $(addprefix $(CXX_OBJDIR)/,$(CXX_OBJS))

#--------- Add to list of known file suffixes
.SUFFIXES: .$(cxxSrcSuf) .$(fSrcSuf) .$(c++SrcSuf) .$(cSrcSuf)

.phony: clean tidy

all: $(FORT_OBJS_W_DIR) $(CXX_OBJS_W_DIR) $(PIXIE)

$(FORT_OBJS_W_DIR): | $(FORT_OBJDIR)

$(FORT_OBJDIR):
	mkdir -p $(FORT_OBJDIR)

$(FORT_OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

$(CXX_OBJS_W_DIR): | $(CXX_OBJDIR)

$(CXX_OBJDIR):
	mkdir -p $(CXX_OBJDIR)

$(CXX_OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

#----------- link all created objects together
#----------- to create pixie_ldf_c program
$(PIXIE): $(FORT_OBJS_W_DIR) $(CXX_OBJS_W_DIR) $(LIBS)
	$(LINK.o) $^ -o $@ $(LDLIBS)

Test: $(FORT_OBJS_W_DIR) $(CXX_OBJS_W_DIR) $(LIBS)
	$(LINK.o) $^ -o $@ $(LDLIBS)

#----------- remove all objects, core and .so file
clean:
	@echo "Cleaning up..."
	@rm -f $(CXX_OBJDIR)/*.o $(PIXIE) core *~ src/*~ include/*~ scan/*~
#$(FORT_OBJDIR)/*.o 
tidy:
	@echo "Tidying up..."
	@rm -f core *~ src/*~ include/*~ scan/*~
