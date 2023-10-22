add_executable(test_operators test_operators.cpp)
target_link_libraries(test_operators PRIVATE tinned ${SYMENGINE_LIBRARIES})
add_test(NAME test_operators COMMAND test_operators)

