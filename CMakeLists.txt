cmake_minimum_required(VERSION 3.6)
project(EllipticCurveMethod)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -o -pthread -lgmp -lm")


set(SOURCE_FILES main.c dataStructures.h ECarithmetics.c)
add_executable(ECM ${SOURCE_FILES})