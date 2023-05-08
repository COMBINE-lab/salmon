# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


set(CPACK_BUILD_SOURCE_DIRS "/home/zzare/salmon;/home/zzare/salmon/build")
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
set(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_FILE "/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_SUMMARY "Salmon built using CMake")
set(CPACK_GENERATOR "TGZ")
set(CPACK_INSTALL_CMAKE_PROJECTS "/home/zzare/salmon/build;Salmon;ALL;/")
set(CPACK_INSTALL_PREFIX "/home/zzare/salmon")
set(CPACK_MODULE_PATH "/home/zzare/salmon/cmake/Modules/")
set(CPACK_NSIS_DISPLAY_NAME "Salmon-1.9.0 1.9.0")
set(CPACK_NSIS_INSTALLER_ICON_CODE "")
set(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
set(CPACK_NSIS_PACKAGE_NAME "Salmon-1.9.0 1.9.0")
set(CPACK_NSIS_UNINSTALL_NAME "Uninstall")
set(CPACK_OBJCOPY_EXECUTABLE "/usr/bin/objcopy")
set(CPACK_OBJDUMP_EXECUTABLE "/usr/bin/objdump")
set(CPACK_OUTPUT_CONFIG_FILE "/home/zzare/salmon/build/CPackConfig.cmake")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
set(CPACK_PACKAGE_DESCRIPTION_FILE "/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Salmon - Wicked-fast RNA-seq isoform quantification using selective alignment")
set(CPACK_PACKAGE_FILE_NAME "Salmon-1.9.0-1.9.0-Linux")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "Salmon-1.9.0 1.9.0")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "Salmon-1.9.0 1.9.0")
set(CPACK_PACKAGE_NAME "Salmon-1.9.0")
set(CPACK_PACKAGE_RELOCATABLE "true")
set(CPACK_PACKAGE_VENDOR "University of Maryland")
set(CPACK_PACKAGE_VERSION "1.9.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "9")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_READELF_EXECUTABLE "/usr/bin/readelf")
set(CPACK_RESOURCE_FILE_LICENSE "/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/share/cmake/Templates/CPack.GenericLicense.txt")
set(CPACK_RESOURCE_FILE_README "/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/share/cmake/Templates/CPack.GenericDescription.txt")
set(CPACK_RESOURCE_FILE_WELCOME "/home/linuxbrew/.linuxbrew/Cellar/cmake/3.25.1/share/cmake/Templates/CPack.GenericWelcome.txt")
set(CPACK_SET_DESTDIR "OFF")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/src/PCA.cpp;/src/PCAUtils.cpp;/build/;/scripts/AggregateToGeneLevel.py;/scripts/ExpressionTools.py;/scripts/GenerateExpressionFiles.sh;/scripts/ParseSoftFile.py;/scripts/PlotCorrelation.py;/scripts/junk;/scripts/sfstrace.log;/scripts/SFPipeline.py;/bin/;/lib/;/sample_data/;PublishREADMEToWebsite.sh;/external/;/src/obsolete/;/include/obsolete/;WebsiteHeader.txt;/experimental_configs/;.git/")
set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/zzare/salmon/build/CPackSourceConfig.cmake")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "Salmon-1.9.0-Source")
set(CPACK_SYSTEM_NAME "Linux")
set(CPACK_THREADS "1")
set(CPACK_TOPLEVEL_TAG "Linux")
set(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/home/zzare/salmon/build/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()
