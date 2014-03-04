################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/code/src/seis/GR.c \
../code/code/src/seis/Helmstetter.c \
../code/code/src/seis/WellsCoppersmith.c \
../code/code/src/seis/background_rate.c \
../code/code/src/seis/cmbopt.c \
../code/code/src/seis/decluster.c \
../code/code/src/seis/soumod1.c 

OBJS += \
./code/code/src/seis/GR.o \
./code/code/src/seis/Helmstetter.o \
./code/code/src/seis/WellsCoppersmith.o \
./code/code/src/seis/background_rate.o \
./code/code/src/seis/cmbopt.o \
./code/code/src/seis/decluster.o \
./code/code/src/seis/soumod1.o 

C_DEPS += \
./code/code/src/seis/GR.d \
./code/code/src/seis/Helmstetter.d \
./code/code/src/seis/WellsCoppersmith.d \
./code/code/src/seis/background_rate.d \
./code/code/src/seis/cmbopt.d \
./code/code/src/seis/decluster.d \
./code/code/src/seis/soumod1.d 


# Each subdirectory must supply rules for building sources it contributes
code/code/src/seis/%.o: ../code/code/src/seis/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -pg -Wall -c -fmessage-length=0 -std=c99 -fopenmp -pg -lxml2 -I/usr/include/libxml2 -lxslt -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


