#######################################################
# trainVEGASneuralNet
#######################################################

add_executable(sst_gate_test sst_gate_test.cpp)

target_link_libraries(sst_gate_test
		${ROOT_LIBRARIES}
		${ROBAST_LIBRARY}
		)

#target_link_libraries(trainVEGASneuralNet trainVEGASneuralNetLib)

install(TARGETS sst_gate_test DESTINATION bin)