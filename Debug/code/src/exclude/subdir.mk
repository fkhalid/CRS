################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/src/exclude/forecast_stepG.c 

OBJS += \
./code/src/exclude/forecast_stepG.o 

C_DEPS += \
./code/src/exclude/forecast_stepG.d 


# Each subdirectory must supply rules for building sources it contributes
code/src/exclude/%.o: ../code/src/exclude/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -pg -Wall -c -fmessage-length=0 -std=c99 -fopenmp -pg -lxml2 -I/usr/include/libxml2 -lxslt -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


