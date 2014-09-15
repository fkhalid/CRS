################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/src/figures_CompPaper.c 

OBJS += \
./code/src/figures_CompPaper.o 

C_DEPS += \
./code/src/figures_CompPaper.d 


# Each subdirectory must supply rules for building sources it contributes
code/src/%.o: ../code/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


