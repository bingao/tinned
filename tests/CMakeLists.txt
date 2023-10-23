add_executable(test_operators test_operators.cpp)
target_link_libraries(test_operators PRIVATE tinned ${SYMENGINE_LIBRARIES})
add_test(NAME test_operators COMMAND test_operators)

add_executable(test_visitors test_visitors.cpp)
target_link_libraries(test_visitors PRIVATE tinned ${SYMENGINE_LIBRARIES})
add_test(NAME test_visitors COMMAND test_visitors)
