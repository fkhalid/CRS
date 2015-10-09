################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
other/bccbug.c \
other/complex.c \
other/nrutil.c 

OBJS += \
./nr_other/bccbug.o \
./nr_other/complex.o \
./nr_other/nrutil.o 

C_DEPS += \
./nr_other/bccbug.d \
./nr_other/complex.d \
./nr_other/nrutil.d 


# Each subdirectory must supply rules for building sources it contributes
nr_other/bccbug.o: other/bccbug.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -D_no_numericalrecipes -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

nr_other/complex.o: other/complex.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -D_no_numericalrecipes -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

nr_other/nrutil.o: other/nrutil.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -D_no_numericalrecipes -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


