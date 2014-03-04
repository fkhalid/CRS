################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/code/src/geom/convert_geometry.c \
../code/code/src/geom/coord_trafos.c \
../code/code/src/geom/dft2d.c \
../code/code/src/geom/dist2fault.c \
../code/code/src/geom/find_gridpoints.c 

OBJS += \
./code/code/src/geom/convert_geometry.o \
./code/code/src/geom/coord_trafos.o \
./code/code/src/geom/dft2d.o \
./code/code/src/geom/dist2fault.o \
./code/code/src/geom/find_gridpoints.o 

C_DEPS += \
./code/code/src/geom/convert_geometry.d \
./code/code/src/geom/coord_trafos.d \
./code/code/src/geom/dft2d.d \
./code/code/src/geom/dist2fault.d \
./code/code/src/geom/find_gridpoints.d 


# Each subdirectory must supply rules for building sources it contributes
code/code/src/geom/%.o: ../code/code/src/geom/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -pg -Wall -c -fmessage-length=0 -std=c99 -fopenmp -pg -lxml2 -I/usr/include/libxml2 -lxslt -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


