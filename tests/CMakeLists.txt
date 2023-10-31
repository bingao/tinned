# Adds Catch2::Catch2, from https://gitlab.com/CLIUtils/modern-cmake.git
FetchContent_Declare(
    catch
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v2.13.6)
FetchContent_MakeAvailable(catch)

add_executable(test_operators test_operators.cpp)
target_link_libraries(test_operators PRIVATE tinned ${SYMENGINE_LIBRARIES} Catch2::Catch2)
add_test(NAME test_operators COMMAND test_operators)

#add_executable(test_visitors test_visitors.cpp)
#target_link_libraries(test_visitors PRIVATE tinned ${SYMENGINE_LIBRARIES} Catch2::Catch2)
#add_test(NAME test_visitors COMMAND test_visitors)
