set(Boost_USE_MULTITHREADED ON)

enable_testing()

include_directories(
   ${Boost_INCLUDE_DIRS}
)
add_executable(tests
   test_entry.cpp
   test_signal.cpp
   test_fourier.cpp
   test_convolve.cpp
)
target_include_directories(tests PRIVATE
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>/tempmeta
)
target_link_libraries(tests
   Eigen3::Eigen
   ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

target_compile_features(tests PRIVATE cxx_std_17)
add_test(NAME dsp_tests COMMAND tests -o report_tests.xml -r junit)

