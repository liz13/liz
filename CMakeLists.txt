cmake_minimum_required(VERSION 3.0.0)

project(CoupledDD VERSION 0.0.1 LANGUAGES CXX)

include(cmake/bob.cmake)

bob_begin_package()

set(Scorec_PREFIX "" CACHE PATH "Scorec installation directory")

bob_begin_cxx_flags()
bob_cxx11_flags()
bob_end_cxx_flags()

include(cmake/dependencies.cmake)
add_subdirectory(src)

bob_end_package()
