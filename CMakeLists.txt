# use cmake 2.8 or later
cmake_minimum_required (VERSION 3.2)

project(minostest)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
if (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
endif()

set(CMAKE_MODULE_PATH $ENV{CANVAS_DIR}/Modules
                      $ENV{CETBUILDTOOLS_DIR}/Modules
                      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

set_install_root()

cet_set_compiler_flags(DIAGS CAUTIOUS
  #WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -Wno-pedantic -Wno-unused-local-typedefs -Wno-narrowing -Wno-deprecated-declarations -Wno-unused-variable -Wno-missing-braces -Wno-undefined-var-template
  EXTRA_CXX_FLAGS
)

cet_enable_asserts()

# issue 5589: if we are not in a debugging build, disable the dispatch of debug messages in messagefacility
string(TOUPPER ${CMAKE_BUILD_TYPE} BTYPE_UC )
if( ${BTYPE_UC} STREQUAL "OPT" OR
    ${BTYPE_UC} STREQUAL "PROF" OR
    ${BTYPE_UC} STREQUAL "RELEASE" OR
    ${BTYPE_UC} STREQUAL "MINSIZEREL" )
  message ( STATUS "Disabling LOG_DEBUG from messagefacility library" )
  add_definitions(-DML_NDEBUG)
endif()

add_definitions( -DNOVACMAKE )

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
# put products in this list that you want to have set up with
# their environment by ups
find_ups_product( cetbuildtools  v4_07_03  )
find_ups_product( cafanacore     v01.00    )
find_ups_product( castxml        v0_00_00_f1 )
find_ups_product( dk2nudata      v01_06_00 )
find_ups_product( dk2nugenie     v01_06_00 )
find_ups_product( eigen          v3.3.5    )
find_ups_product( gcc            v4_9_3    )
find_ups_product( genie          v2_10_00  )
find_ups_product( gsl            v1_16     )
find_ups_product( libwda         v2_20_1   )
find_ups_product( novadaq        v10_00_00 )
find_ups_product( novarwgt       v00.13    )
find_ups_product( nuevdb         v1_00_00  )
find_ups_product( nutools        v1_21_00  )
find_ups_product( nugen          v1_02_00  )
find_ups_product( nug4           v1_02_00  )
find_ups_product( opencv         v3_4_2    )
find_ups_product( ppfx           v01_09    )
find_ups_product( pygccxml       v1_9_1    )
find_ups_product( osclib         v00.03    )
find_ups_product( protobuf       v3_3_1a   )
find_ups_product( stan           v2.18.0   )
find_ups_product( stan_math      v2.18.0   )
find_ups_product( tensorflow     v1_3_0b   )


find_ups_root( root v5_00_00 )

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)




# source

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make(LIBRARY_NAME     MinosTest
         LIB_LIBRARIES    ${ROOT_BASIC_LIB_LIST}
	 MODULE_LIBRARIES ${CETLIB} ${CETLIB_EXCEPT}
	 		  ${ROOT_BASIC_LIB_LIST}
          BASENAME_ONLY
         )
install_headers()
install_source()

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)

