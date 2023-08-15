cmake_minimum_required(VERSION 3.24)
project(compositeILM)

set(CMAKE_CXX_STANDARD 17)
#find_package(Armadillo REQUIRED)

# Find Armadillo package
find_package(Armadillo REQUIRED)
# Add your source files
add_executable(compositeILM main.cpp ILM.cpp MCMC.cpp ILM.h MCMC.h generate_individuals.cpp generate_individuals.h simulate.cpp simulate.h Inference.cpp Inference.cpp Inference.h Clustering.cpp Clustering.cpp Clustering.h HamiltonianMC.cpp HamiltonianMC.cpp HamiltonianMC.h Inference_composite.cpp Inference_composite.cpp Inference_composite.h MCMC_composite.cpp MCMC_composite.h)
# Create the executable and link libraries
#add_executable(compositeILM ${SOURCE_FILES})
target_link_libraries(compositeILM PRIVATE ${ARMADILLO_LIBRARIES} /usr/local/lib/libcppad_lib.dylib)
target_include_directories(compositeILM PRIVATE ${ARMADILLO_INCLUDE_DIRS})
include_directories(/usr/local/include)
#target_link_libraries(compositeILM /usr/local/lib/libcppad_lib.dylib)