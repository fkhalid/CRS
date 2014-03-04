################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/code/src/inp_out/print_output.c \
../code/code/src/inp_out/propagate_results.c \
../code/code/src/inp_out/read_crust.c \
../code/code/src/inp_out/read_csep_template.c \
../code/code/src/inp_out/read_eqkfm.c \
../code/code/src/inp_out/read_focmec.c \
../code/code/src/inp_out/read_inputfile.c \
../code/code/src/inp_out/read_matrix.c \
../code/code/src/inp_out/read_param.c \
../code/code/src/inp_out/read_zmap.c \
../code/code/src/inp_out/write_csep_forecast.c 

OBJS += \
./code/code/src/inp_out/print_output.o \
./code/code/src/inp_out/propagate_results.o \
./code/code/src/inp_out/read_crust.o \
./code/code/src/inp_out/read_csep_template.o \
./code/code/src/inp_out/read_eqkfm.o \
./code/code/src/inp_out/read_focmec.o \
./code/code/src/inp_out/read_inputfile.o \
./code/code/src/inp_out/read_matrix.o \
./code/code/src/inp_out/read_param.o \
./code/code/src/inp_out/read_zmap.o \
./code/code/src/inp_out/write_csep_forecast.o 

C_DEPS += \
./code/code/src/inp_out/print_output.d \
./code/code/src/inp_out/propagate_results.d \
./code/code/src/inp_out/read_crust.d \
./code/code/src/inp_out/read_csep_template.d \
./code/code/src/inp_out/read_eqkfm.d \
./code/code/src/inp_out/read_focmec.d \
./code/code/src/inp_out/read_inputfile.d \
./code/code/src/inp_out/read_matrix.d \
./code/code/src/inp_out/read_param.d \
./code/code/src/inp_out/read_zmap.d \
./code/code/src/inp_out/write_csep_forecast.d 


# Each subdirectory must supply rules for building sources it contributes
code/code/src/inp_out/%.o: ../code/code/src/inp_out/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -pg -Wall -c -fmessage-length=0 -std=c99 -fopenmp -pg -lxml2 -I/usr/include/libxml2 -lxslt -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


