################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../AlignedHelix.cpp \
../Alignment.cpp \
../AlignmentGenerator.cpp \
../BasePair.cpp \
../CompetingHelix.cpp \
../EvolModel.cpp \
../Helix.cpp \
../HelixCore.cpp \
../HelixFinder.cpp \
../HelixGroup.cpp \
../InterestingRegion.cpp \
../SeqHelix.cpp \
../ShuffledAlignment.cpp \
../StatsWrapper.cpp \
../TransatMain.cpp \
../Tree.cpp \
../UTMatrix.cpp \
../Utilities.cpp 

OBJS += \
./AlignedHelix.o \
./Alignment.o \
./AlignmentGenerator.o \
./BasePair.o \
./CompetingHelix.o \
./EvolModel.o \
./Helix.o \
./HelixCore.o \
./HelixFinder.o \
./HelixGroup.o \
./InterestingRegion.o \
./SeqHelix.o \
./ShuffledAlignment.o \
./StatsWrapper.o \
./TransatMain.o \
./Tree.o \
./UTMatrix.o \
./Utilities.o 

CPP_DEPS += \
./AlignedHelix.d \
./Alignment.d \
./AlignmentGenerator.d \
./BasePair.d \
./CompetingHelix.d \
./EvolModel.d \
./Helix.d \
./HelixCore.d \
./HelixFinder.d \
./HelixGroup.d \
./InterestingRegion.d \
./SeqHelix.d \
./ShuffledAlignment.d \
./StatsWrapper.d \
./TransatMain.d \
./Tree.d \
./UTMatrix.d \
./Utilities.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

ShuffledAlignment.o: ../ShuffledAlignment.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DRNAZ_SHUFFLER_LOC=\"$(RNAZ_SHUFFLER_LOC)\" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"ShuffledAlignment.d" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


