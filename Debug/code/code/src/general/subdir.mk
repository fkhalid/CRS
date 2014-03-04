################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/code/src/general/CRS_LogLikelihood.c \
../code/code/src/general/calculateDCFSperturbed.c \
../code/code/src/general/eqkfm_copy.c \
../code/code/src/general/find_timesteps.c \
../code/code/src/general/forecast_stepG.c \
../code/code/src/general/mem_mgmt.c \
../code/code/src/general/setup.c \
../code/code/src/general/struct_conversions.c 

OBJS += \
./code/code/src/general/CRS_LogLikelihood.o \
./code/code/src/general/calculateDCFSperturbed.o \
./code/code/src/general/eqkfm_copy.o \
./code/code/src/general/find_timesteps.o \
./code/code/src/general/forecast_stepG.o \
./code/code/src/general/mem_mgmt.o \
./code/code/src/general/setup.o \
./code/code/src/general/struct_conversions.o 

C_DEPS += \
./code/code/src/general/CRS_LogLikelihood.d \
./code/code/src/general/calculateDCFSperturbed.d \
./code/code/src/general/eqkfm_copy.d \
./code/code/src/general/find_timesteps.d \
./code/code/src/general/forecast_stepG.d \
./code/code/src/general/mem_mgmt.d \
./code/code/src/general/setup.d \
./code/code/src/general/struct_conversions.d 


# Each subdirectory must supply rules for building sources it contributes
code/code/src/general/%.o: ../code/code/src/general/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -pg -Wall -c -fmessage-length=0 -std=c99 -fopenmp -pg -lxml2 -I/usr/include/libxml2 -lxslt -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


