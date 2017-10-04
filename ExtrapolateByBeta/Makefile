# Build the extrap by limit finder

ROOTCFLAGS	= $(shell root-config --cflags)
ROOTLIBS	= $(shell root-config --libs)
ROOTGLIBS	= $(shell root-config --glibs)

CXX		= gcc
CXXFLAGS	=-I$(ROOTSYS)/include -O -Wall -fPIC
LD		= gcc
LDFLAGS		= -g
SOFLAGS		= -shared

COMMONLIM	= ../LimitCommonCode
COMMONUTILS	= ../CommonLimitUtils
WILD		= ..
CXXFLAGS	+= $(ROOTCFLAGS) -I$(COMMONLIM) -I$(COMMONUTILS) -I$(WILD)
LIBS    = $(ROOTLIBS) -L$(TMVASYS)/lib -lTMVA $(shell root-config --libs) -lMLP -lXMLIO -lTreePlayer -lstdc++ -lRooFitCore -lRooStats -lRooFit -lMinuit -lFoam -lMathMore
GLIBS		= $(ROOTGLIBS)

OBJS		= main.o Lxy_weight_calculator.o muon_tree_processor.o limitSetting.o run_ABCD.o HypoTestInvTool.o

ExtrapolateByBeta:	$(OBJS)
	$(CXX) -o $@ $(OBJS) $(CXXFLAGS) $(LIBS)

main.o : main.cxx $(COMMONLIM)/extrap_file_wrapper.h $(COMMONLIM)/limitSetting.h
	$(CXX) -c main.cxx $(CXXFLAGS)

Lxy_weight_calculator.o : Lxy_weight_calculator.cxx $(COMMONLIM)/extrap_file_wrapper.h $(COMMONLIM)/limitSetting.h
	$(CXX) -c Lxy_weight_calculator.cxx $(CXXFLAGS)

muon_tree_processor.o : muon_tree_processor.cxx $(COMMONLIM)/extrap_file_wrapper.h $(COMMONLIM)/limitSetting.h
	$(CXX) -c muon_tree_processor.cxx $(CXXFLAGS)

limitSetting.o : $(COMMONLIM)/limitSetting.cxx $(COMMONLIM)/limitSetting.h $(COMMONLIM)/extrap_file_wrapper.h $(COMMONLIM)/limit_output_file.h $(COMMONLIM)/SimulABCD.h $(COMMONLIM)/limit_datastructures.h
	$(CXX) -c $(COMMONLIM)/limitSetting.cxx $(CXXFLAGS)
	
run_ABCD.o : $(COMMONLIM)/run_ABCD.cxx $(COMMONLIM)/HypoTestInvTool.h $(COMMONLIM)/SimulABCD.h
	$(CXX) -c $(COMMONLIM)/run_ABCD.cxx $(CXXFLAGS)

HypoTestInvTool.o : $(COMMONLIM)/HypoTestInvTool.cxx $(COMMONLIM)/HypoTestInvTool.h
	$(CXX) -c $(COMMONLIM)/HypoTestInvTool.cxx $(CXXFLAGS)

	# suffix rule
.cxx.o:
	$(CXX) -c $(CXXFLAGS) $< #$<

# suffix rule
.C.o:
	$(CXX) -c $(CXXFLAGS) $< #$<

# clean
clean:
	rm -f *~ *.o *.o~ core run
