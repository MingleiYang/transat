################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/AlignedHelix.cpp \
../src/Alignment.cpp \
../src/AlignmentGenerator.cpp \
../src/BasePair.cpp \
../src/CompetingHelix.cpp \
../src/EvolModel.cpp \
../src/Helix.cpp \
../src/HelixCore.cpp \
../src/HelixFinder.cpp \
../src/HelixGroup.cpp \
../src/InterestingRegion.cpp \
../src/SeqHelix.cpp \
../src/ShuffledAlignment.cpp \
../src/StatsWrapper.cpp \
../src/Tree.cpp \
../src/UTMatrix.cpp \
../src/Utilities.cpp \
../src/get_all_statistics.cpp 

OBJS += \
./src/AlignedHelix.o \
./src/Alignment.o \
./src/AlignmentGenerator.o \
./src/BasePair.o \
./src/CompetingHelix.o \
./src/EvolModel.o \
./src/Helix.o \
./src/HelixCore.o \
./src/HelixFinder.o \
./src/HelixGroup.o \
./src/InterestingRegion.o \
./src/SeqHelix.o \
./src/ShuffledAlignment.o \
./src/StatsWrapper.o \
./src/Tree.o \
./src/UTMatrix.o \
./src/Utilities.o \
./src/get_all_statistics.o 

CPP_DEPS += \
./src/AlignedHelix.d \
./src/Alignment.d \
./src/AlignmentGenerator.d \
./src/BasePair.d \
./src/CompetingHelix.d \
./src/EvolModel.d \
./src/Helix.d \
./src/HelixCore.d \
./src/HelixFinder.d \
./src/HelixGroup.d \
./src/InterestingRegion.d \
./src/SeqHelix.d \
./src/ShuffledAlignment.d \
./src/StatsWrapper.d \
./src/Tree.d \
./src/UTMatrix.d \
./src/Utilities.d \
./src/get_all_statistics.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


