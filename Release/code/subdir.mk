################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/main.c 

OBJS += \
./code/main.o 

C_DEPS += \
./code/main.d 


# Each subdirectory must supply rules for building sources it contributes
code/%.o: ../code/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I/home/des/camcat/masterdes_sec21/NumericalRecipes/2.11/C_211/other -I/home/des/camcat/masterdes_sec21/NumericalRecipes/2.11/C_211/recipes -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


