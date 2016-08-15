##
# Print some post install messages for the user
##
message("\n\n")
message("Installation complete. Please ensure the following paths are set properly.")	
message("==========================================================================")
#if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#	message("Fixing library names with install_name_tool")
#	execute_process(COMMAND install_name_tool -add_rpath ${CMAKE_INSTALL_PREFIX}/bin ${CMAKE_INSTALL_PREFIX}/bin/salmon)
#	execute_process(COMMAND install_name_tool -add_rpath ${CMAKE_INSTALL_PREFIX}/lib ${CMAKE_INSTALL_PREFIX}/bin/salmon)
#	execute_process(COMMAND install_name_tool -add_rpath @executable_path ${CMAKE_INSTALL_PREFIX}/bin/salmon) 
#endif()
message("Please add ${CMAKE_INSTALL_PREFIX}/bin to your PATH")
if ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
	message("Please add ${CMAKE_INSTALL_PREFIX}/lib to your DYLD_FALLBACK_LIBRARY_PATH")
else()
	message("Please add ${CMAKE_INSTALL_PREFIX}/lib to your LD_LIBRARY_PATH")
endif()
message("==========================================================================")
