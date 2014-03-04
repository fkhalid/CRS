################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/code/src/okada/dc3d.c \
../code/code/src/okada/okadaDCFS.c \
../code/code/src/okada/prestress.c \
../code/code/src/okada/pscokada.c 

OBJS += \
./code/code/src/okada/dc3d.o \
./code/code/src/okada/okadaDCFS.o \
./code/code/src/okada/prestress.o \
./code/code/src/okada/pscokada.o 

C_DEPS += \
./code/code/src/okada/dc3d.d \
./code/code/src/okada/okadaDCFS.d \
./code/code/src/okada/prestress.d \
./code/code/src/okada/pscokada.d 


# Each subdirectory must supply rules for building sources it contributes
code/code/src/okada/%.o: ../code/code/src/okada/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -pg -Wall -c -fmessage-length=0 -std=c99 -fopenmp -pg -lxml2 -I/usr/include/libxml2 -lxslt -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


