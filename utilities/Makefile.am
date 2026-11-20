AM_LDFLAGS += -L$(top_builddir)/Hadrons
AM_CXXFLAGS += -I$(top_srcdir)

bin_PROGRAMS = \
  HadronsXmlRun              \
  HadronsXmlValidate

EXTRA_PROGRAMS = \
  HadronsContractor          \
  ContractorDense            \
  K2pipiContractor           \
  BatchDeflationBenchmark    \
  HadronsContractorBenchmark \
  HadronsXmlRun              \
  HadronsXmlValidate         \
  HadronsFermionEP64To32     \
  BubbleContractor           \
  A2AMesonFieldSelectHits    \
  A2AMesonFieldCast          \
  FlexibleContractor         \
  A2APixPi                   \
  FlexibleContractorDense    \
  BubbleContractorDense      \
  K2pipiContractorDense      \
  A2AMesonFieldProduct       \
  FlexibleContractorDWC      \
  BubbleContractorDWC        \
  K2pipiContractorDWC        \
  HadronsFermionEP64To32

CLEANFILES = $(EXTRA_PROGRAMS)

BatchDeflationBenchmark_SOURCES = BatchDeflationBenchmark.cpp
BatchDeflationBenchmark_LDADD = -lHadrons -lGrid
BatchDeflationBenchmark_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a

HadronsXmlRun_SOURCES = HadronsXmlRun.cpp
HadronsXmlRun_LDADD   = -lHadrons -lGrid
HadronsXmlRun_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a

HadronsXmlValidate_SOURCES = HadronsXmlValidate.cpp
HadronsXmlValidate_LDADD   = -lHadrons -lGrid
HadronsXmlValidate_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a

HadronsFermionEP64To32_SOURCES  = EigenPackCast.cpp
HadronsFermionEP64To32_CXXFLAGS = $(AM_CXXFLAGS) -DFIN=WilsonImplD::FermionField -DFOUT=WilsonImplF::FermionField
HadronsFermionEP64To32_LDADD    = -lHadrons -lGrid
HadronsFermionEP64To32_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a

HadronsContractor_SOURCES = Contractor.cpp
HadronsContractor_LDADD   = -lHadrons -lGrid
HadronsContractor_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a

ContractorDense_SOURCES = ContractorDense.cpp
ContractorDense_LDADD   = -lHadrons -lGrid

FlexibleContractor_SOURCES = FlexibleContractor.cpp
FlexibleContractor_LDADD   = -lHadrons -lGrid

K2pipiContractor_SOURCES = K2pipiContractor.cpp
K2pipiContractor_LDADD   = -lHadrons -lGrid

BubbleContractor_SOURCES = BubbleContractor.cpp
BubbleContractor_LDADD   = -lHadrons -lGrid

A2AMesonFieldSelectHits_SOURCES = A2AMesonFieldSelectHits.cpp
A2AMesonFieldSelectHits_LDADD   = -lHadrons -lGrid

A2APixPi_SOURCES = A2APixPi.cpp A2ADenseContract.hpp
A2APixPi_LDADD   = -lHadrons -lGrid

FlexibleContractorDense_SOURCES = FlexibleContractorDense.cpp A2ADenseContract.hpp
FlexibleContractorDense_LDADD   = -lHadrons -lGrid

BubbleContractorDense_SOURCES = BubbleContractorDense.cpp A2ADenseContract.hpp
BubbleContractorDense_LDADD   = -lHadrons -lGrid

K2pipiContractorDense_SOURCES = K2pipiContractorDense.cpp A2ADenseContract.hpp
K2pipiContractorDense_LDADD   = -lHadrons -lGrid

A2AMesonFieldProduct_SOURCES = A2AMesonFieldProduct.cpp
A2AMesonFieldProduct_LDADD   = -lHadrons -lGrid

FlexibleContractorDWC_SOURCES = FlexibleContractorDWC.cpp
FlexibleContractorDWC_LDADD   = -lHadrons -lGrid

BubbleContractorDWC_SOURCES = BubbleContractorDWC.cpp
BubbleContractorDWC_LDADD   = -lHadrons -lGrid

K2pipiContractorDWC_SOURCES = K2pipiContractorDWC.cpp
K2pipiContractorDWC_LDADD   = -lHadrons -lGrid

A2AMesonFieldCast_SOURCES = A2AMesonFieldCast.cpp
A2AMesonFieldCast_LDADD   = -lHadrons -lGrid

HadronsContractorBenchmark_SOURCES = ContractorBenchmark.cpp
HadronsContractorBenchmark_LDADD   = -lHadrons -lGrid
HadronsContractorBenchmark_DEPENDENCIES = $(top_builddir)/Hadrons/libHadrons.a
