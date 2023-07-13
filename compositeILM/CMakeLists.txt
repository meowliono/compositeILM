cmake_minimum_required(VERSION 3.24)
project(compositeILM)

set(CMAKE_CXX_STANDARD 17)

add_executable(compositeILM main.cpp ILM.cpp MCMC.cpp ILM.h MCMC.h generate_individuals.cpp generate_individuals.h simulate.cpp simulate.h Inference.cpp Inference.cpp Inference.h)
