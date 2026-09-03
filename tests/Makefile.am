AM_LDFLAGS  += -L$(top_builddir)/Hadrons
AM_CXXFLAGS += -I$(top_srcdir)
 
EXTRA_PROGRAMS = \
  Test_QED_Local            \
  Test_QED                  \
  Test_database             \
  Test_database_concurrency \
  Test_diskvector           \
  Test_em_field             \
  Test_exact_distil         \
  Test_fermion_flow			\
  Test_field_io             \
  Test_free_prop            \
  Test_gradient_flow		\
  Test_hadrons_meson_3pt    \
  Test_hadrons_spectrum     \
  Test_highfreq_stat        \
  Test_result_bundling      \
  Test_sigma_to_nucleon     \
  Test_stoch_distil         \
  Test_xi_to_sigma          \
  Test_emf_gpu              \
  Test_emf_compare          \
  Test_cmof_gpu             \
  Test_cmof_compare         \
  Test_mf_gpu               \
  Test_mf_compare           \
  Test_writer_spread        \
  Test_load_combined_regression \
  Test_mf_multihit          \
  Test_savebinned_gpu       \
  Test_a2a_evalfpe_gpu      \
  Test_coarsedeflation_gpu

CLEANFILES = $(EXTRA_PROGRAMS)

tests-local: $(EXTRA_PROGRAMS)

Test_QED_Local_SOURCES=Test_QED_Local.cpp
Test_QED_Local_LDADD=-lHadrons -lGrid

Test_QED_SOURCES=Test_QED.cpp
Test_QED_LDADD=-lHadrons -lGrid

Test_diskvector_SOURCES=Test_diskvector.cpp
Test_diskvector_LDADD=-lHadrons -lGrid

Test_database_SOURCES=Test_database.cpp
Test_database_LDADD=-lHadrons -lGrid

Test_database_concurrency_SOURCES=Test_database_concurrency.cpp
Test_database_concurrency_LDADD=-lHadrons -lGrid

Test_em_field_SOURCES=Test_em_field.cpp
Test_em_field_LDADD=-lHadrons -lGrid

Test_exact_distil_SOURCES=Test_exact_distil.cpp
Test_exact_distil_LDADD=-lHadrons -lGrid

Test_fermion_flow_SOURCES=Test_fermion_flow.cpp
Test_fermion_flow_LDADD=-lHadrons -lGrid

Test_field_io_SOURCES=Test_field_io.cpp
Test_field_io_LDADD=-lHadrons -lGrid

Test_free_prop_SOURCES=Test_free_prop.cpp
Test_free_prop_LDADD=-lHadrons -lGrid

Test_gradient_flow_SOURCES=Test_gradient_flow.cpp
Test_gradient_flow_LDADD=-lHadrons -lGrid

Test_hadrons_meson_3pt_SOURCES=Test_hadrons_meson_3pt.cpp
Test_hadrons_meson_3pt_LDADD=-lHadrons -lGrid

Test_hadrons_spectrum_SOURCES=Test_hadrons_spectrum.cpp
Test_hadrons_spectrum_LDADD=-lHadrons -lGrid

Test_highfreq_stat_SOURCES=Test_highfreq_stat.cpp
Test_highfreq_stat_LDADD=-lHadrons -lGrid

Test_result_bundling_SOURCES=Test_result_bundling.cpp
Test_result_bundling_LDADD=-lHadrons -lGrid

Test_sigma_to_nucleon_SOURCES=Test_sigma_to_nucleon.cpp
Test_sigma_to_nucleon_LDADD=-lHadrons -lGrid

Test_stoch_distil_SOURCES=Test_stoch_distil.cpp
Test_stoch_distil_LDADD=-lHadrons -lGrid

Test_xi_to_sigma_SOURCES=Test_xi_to_sigma.cpp
Test_xi_to_sigma_LDADD=-lHadrons -lGrid

Test_emf_gpu_SOURCES=Test_emf_gpu.cpp
Test_emf_gpu_LDADD=-lHadrons -lGrid

Test_emf_compare_SOURCES=Test_emf_compare.cpp
Test_emf_compare_LDADD=-lHadrons -lGrid

Test_cmof_gpu_SOURCES=Test_cmof_gpu.cpp
Test_cmof_gpu_LDADD=-lHadrons -lGrid

Test_cmof_compare_SOURCES=Test_cmof_compare.cpp
Test_cmof_compare_LDADD=-lHadrons -lGrid

Test_mf_gpu_SOURCES=Test_mf_gpu.cpp
Test_mf_gpu_LDADD=-lHadrons -lGrid

Test_mf_compare_SOURCES=Test_mf_compare.cpp
Test_mf_compare_LDADD=-lHadrons -lGrid

Test_writer_spread_SOURCES=Test_writer_spread.cpp
Test_writer_spread_LDADD=-lHadrons -lGrid

Test_load_combined_regression_SOURCES=Test_load_combined_regression.cpp
Test_load_combined_regression_LDADD=-lHadrons -lGrid

Test_mf_multihit_SOURCES=Test_mf_multihit.cpp
Test_mf_multihit_LDADD=-lHadrons -lGrid

Test_savebinned_gpu_SOURCES=Test_savebinned_gpu.cpp
Test_savebinned_gpu_LDADD=-lHadrons -lGrid

Test_a2a_evalfpe_gpu_SOURCES=Test_a2a_evalfpe_gpu.cpp
Test_a2a_evalfpe_gpu_LDADD=-lHadrons -lGrid

Test_coarsedeflation_gpu_SOURCES=Test_coarsedeflation_gpu.cpp
Test_coarsedeflation_gpu_LDADD=-lHadrons -lGrid
