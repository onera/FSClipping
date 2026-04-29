# set a variable for the extension of the PYTHONPATH during testing
set(BUILD_PYTHONPATH ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

# thereafter, go over all Python modules one by one
add_subdirectory(FSClipping)
