################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/src/util/eispack.c \
../code/src/util/error.c \
../code/src/util/files.c \
../code/src/util/fit_splines.c \
../code/src/util/gasdev.c \
../code/src/util/gaussj.c \
../code/src/util/hash.c \
../code/src/util/interp_quad.c \
../code/src/util/jacobi.c \
../code/src/util/lnsrch.c \
../code/src/util/merge.c \
../code/src/util/moreutil.c \
../code/src/util/mscorr.c \
../code/src/util/nrutil.c \
../code/src/util/pnpoly.c \
../code/src/util/ran1.c \
../code/src/util/roots3.c \
../code/src/util/spline.c \
../code/src/util/splines_eqkfm.c 

OBJS += \
./code/src/util/eispack.o \
./code/src/util/error.o \
./code/src/util/files.o \
./code/src/util/fit_splines.o \
./code/src/util/gasdev.o \
./code/src/util/gaussj.o \
./code/src/util/hash.o \
./code/src/util/interp_quad.o \
./code/src/util/jacobi.o \
./code/src/util/lnsrch.o \
./code/src/util/merge.o \
./code/src/util/moreutil.o \
./code/src/util/mscorr.o \
./code/src/util/nrutil.o \
./code/src/util/pnpoly.o \
./code/src/util/ran1.o \
./code/src/util/roots3.o \
./code/src/util/spline.o \
./code/src/util/splines_eqkfm.o 

C_DEPS += \
./code/src/util/eispack.d \
./code/src/util/error.d \
./code/src/util/files.d \
./code/src/util/fit_splines.d \
./code/src/util/gasdev.d \
./code/src/util/gaussj.d \
./code/src/util/hash.d \
./code/src/util/interp_quad.d \
./code/src/util/jacobi.d \
./code/src/util/lnsrch.d \
./code/src/util/merge.d \
./code/src/util/moreutil.d \
./code/src/util/mscorr.d \
./code/src/util/nrutil.d \
./code/src/util/pnpoly.d \
./code/src/util/ran1.d \
./code/src/util/roots3.d \
./code/src/util/spline.d \
./code/src/util/splines_eqkfm.d 


# Each subdirectory must supply rules for building sources it contributes
code/src/util/%.o: ../code/src/util/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


