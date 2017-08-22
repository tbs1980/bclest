add_executable(chain_test chain_test.cc)
target_link_libraries(chain_test gtest_main)

add_test(chain_test chain_test)
