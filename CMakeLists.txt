cmake_minimum_required(VERSION 3.14)

include(ExternalProject)

# set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

project(SMG VERSION 1.0.0)
add_executable(smg main.c)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

# minigraph
############
message(STATUS "minigraph will be built from source")
ExternalProject_Add(minigraph
  GIT_REPOSITORY https://github.com/lh3/minigraph.git
  GIT_TAG 07c212dbc2559607772c69152319128c608df373
  BUILD_IN_SOURCE 1
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND "make"
  INSTALL_COMMAND ""
  )
ExternalProject_Get_Property(minigraph SOURCE_DIR)
SET(MG_INCLUDE_DIR ${SOURCE_DIR})

############################################
### final setup (includes and libraries) ###
############################################

add_dependencies(smg minigraph)

target_include_directories(smg
  PRIVATE ${MG_INCLUDE_DIR}
)


target_link_libraries(smg
  PUBLIC ${MG_INCLUDE_DIR}/kalloc.o
  PUBLIC ${MG_INCLUDE_DIR}/kthread.o
  PUBLIC ${MG_INCLUDE_DIR}/algo.o
  PUBLIC ${MG_INCLUDE_DIR}/sys.o
  PUBLIC ${MG_INCLUDE_DIR}/gfa-base.o
  PUBLIC ${MG_INCLUDE_DIR}/gfa-io.o
  PUBLIC ${MG_INCLUDE_DIR}/gfa-aug.o
  PUBLIC ${MG_INCLUDE_DIR}/gfa-bbl.o
  PUBLIC ${MG_INCLUDE_DIR}/gfa-ed.o
  PUBLIC ${MG_INCLUDE_DIR}/sketch.o
  PUBLIC ${MG_INCLUDE_DIR}/misc.o
  PUBLIC ${MG_INCLUDE_DIR}/bseq.o
  PUBLIC ${MG_INCLUDE_DIR}/options.o
  PUBLIC ${MG_INCLUDE_DIR}/shortk.o
  PUBLIC ${MG_INCLUDE_DIR}/miniwfa.o
  PUBLIC ${MG_INCLUDE_DIR}/index.o
  PUBLIC ${MG_INCLUDE_DIR}/lchain.o
  PUBLIC ${MG_INCLUDE_DIR}/gchain1.o
  PUBLIC ${MG_INCLUDE_DIR}/galign.o
  PUBLIC ${MG_INCLUDE_DIR}/gcmisc.o
  PUBLIC ${MG_INCLUDE_DIR}/map-algo.o
  PUBLIC ${MG_INCLUDE_DIR}/cal_cov.o
  PUBLIC ${MG_INCLUDE_DIR}/format.o
  PUBLIC ${MG_INCLUDE_DIR}/gmap.o
  PUBLIC ${MG_INCLUDE_DIR}/ggsimple.o
  PUBLIC ${MG_INCLUDE_DIR}/ggen.o
  PUBLIC ${MG_INCLUDE_DIR}/asm-call.o
  # PUBLIC GCSA
  # PUBLIC SDSL
  PUBLIC z
  PUBLIC m
  PUBLIC pthread
  # PUBLIC OpenMP::OpenMP_CXX
)
