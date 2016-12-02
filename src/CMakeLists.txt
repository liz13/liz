add_executable(driver driver.cpp)
add_executable(driver_edge driver_edge.cpp)
target_link_libraries(driver SCOREC::core)
target_link_libraries(driver_edge SCOREC::core)

bob_end_subdir()
