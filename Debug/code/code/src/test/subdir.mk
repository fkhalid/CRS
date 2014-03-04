################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/code/src/test/d_background_rates.c \
../code/code/src/test/grid_variability.c 

OBJS += \
./code/code/src/test/d_background_rates.o \
./code/code/src/test/grid_variability.o 

C_DEPS += \
./code/code/src/test/d_background_rates.d \
./code/code/src/test/grid_variability.d 


# Each subdirectory must supply rules for building sources it contributes
code/code/src/test/%.o: ../code/code/src/test/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -pg -Wall -c -fmessage-length=0 -std=c99 -fopenmp -pg -lxml2 -I/usr/include/libxml2 -lxslt -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


