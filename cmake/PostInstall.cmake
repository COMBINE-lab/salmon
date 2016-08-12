##
# Print some post install messages for the user
##
message("\n\n")
message("Installation complete. Please ensure the following paths are set properly.")	
message("==========================================================================")
message("Please add ${CMAKE_INSTALL_PREFIX}/bin to your PATH")
if ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
	message("Please add ${CMAKE_INSTALL_PREFIX}/lib to your DYLD_FALLBACK_LIBRARY_PATH")
else()
	message("Please add ${CMAKE_INSTALL_PREFIX}/lib to your LD_LIBRARY_PATH")
endif()
message("==========================================================================")
