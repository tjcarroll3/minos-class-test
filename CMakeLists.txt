
IF( $ENV{MINOSDST} )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# art_make does not know how to handle _plugin.cc files yet, so exclude
# them from being built there.  Make a list of all those files to avoid
# having a huge list in the art_make call
set(EXCLUDED_PLUGINS CMFDSTToEventList_module.cc
   )

art_make(LIBRARY_NAME     CovarianceMatrixFitMinosmodules
         LIB_LIBRARIES    CovarianceMatrixFitcore
                          CovarianceMatrixFitutilities
                          CovarianceMatrixFitdataProducts
                          Utilities
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${MF_MESSAGELOGGER}
                          ${MF_UTILITIES}
                          ${ROOT_BASIC_LIB_LIST}
                          ${FHICLCPP}
         MODULE_LIBRARIES CovarianceMatrixFitdataProducts
                          CovarianceMatrixFitcore
                          CovarianceMatrixFitutilities
                          CovarianceMatrixFitmodules
                          Metadata
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_PERSISTENCY_COMMON}
                          ${MF_MESSAGELOGGER}
                          ${MF_UTILITIES}
                          ${CETLIB} ${CETLIB_EXCEPT}
                          ${ROOT_BASIC_LIB_LIST}
                          ${FHICLCPP}
          EXCLUDE         ${EXCLUDED_PLUGINS}
          BASENAME_ONLY
         )

cet_find_library( BOOST_FILESYSTEM  NAMES boost_filesystem  PATHS ENV BOOST_LIB NO_DEFAULT_PATH )



#  simple_plugin(CMFDSTToEventList "module"
#                ${CAFANA_BASIC_LIB_LIST}
#                3FlavorAnaCuts
#                3FlavorAnaSysts
#                3FlavorAnaVars
#                NuXAnaCuts
#                NuXAnaVars
#                NuXAnaSysts
#                StandardRecord
#                StandardRecordProxy
#                CovarianceMatrixFitcore
#                CovarianceMatrixFitutilities
#                CovarianceMatrixFitdataProducts
#                ${ART_FRAMEWORK_PRINCIPAL}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#                ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#                ${ART_FRAMEWORK_SERVICES_REGISTRY}
#                ${ART_PERSISTENCY_COMMON}
#                ${MF_MESSAGELOGGER}
#                ${MF_UTILITIES}
#                ${CETLIB} ${CETLIB_EXCEPT}
#                ${FHICLCPP}
#                BASENAME_ONLY
#               )



install_source()

ENDIF()
