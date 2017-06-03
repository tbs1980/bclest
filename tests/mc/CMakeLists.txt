add_executable(chain_test chain_test.cc)
target_link_libraries(chain_test ${GTEST_BOTH_LIBRARIES})

add_test(chain_test chain_test)
