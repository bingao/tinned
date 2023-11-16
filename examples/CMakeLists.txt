add_executable(eval_exch_corr eval_exch_corr.cpp)
target_link_libraries(eval_exch_corr PRIVATE tinned ${SYMENGINE_LIBRARIES})
