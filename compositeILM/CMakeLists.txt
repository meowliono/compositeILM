cmake_minimum_required(VERSION 3.24)
project(compositeILM)

set(CMAKE_CXX_STANDARD 17)
#find_package(Armadillo REQUIRED)
include_directories(/usr/local/include)
add_executable(compositeILM main.cpp ILM.cpp MCMC.cpp ILM.h MCMC.h generate_individuals.cpp generate_individuals.h simulate.cpp simulate.h Inference.cpp Inference.cpp Inference.h Clustering.cpp Clustering.cpp Clustering.h HamiltonianMC.cpp HamiltonianMC.cpp HamiltonianMC.h Inference_composite.cpp Inference_composite.cpp Inference_composite.h MCMC_composite.cpp MCMC_composite.h)
target_link_libraries(compositeILM /usr/local/lib/libcppad_lib.dylib)