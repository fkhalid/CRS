################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/src/test/d_background_rates.c \
../code/src/test/grid_variability.c \
../code/src/test/tests.c 

OBJS += \
./code/src/test/d_background_rates.o \
./code/src/test/grid_variability.o \
./code/src/test/tests.o 

C_DEPS += \
./code/src/test/d_background_rates.d \
./code/src/test/grid_variability.d \
./code/src/test/tests.d 


# Each subdirectory must supply rules for building sources it contributes
code/src/test/%.o: ../code/src/test/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	mpicc -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


