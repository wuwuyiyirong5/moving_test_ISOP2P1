################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../ISOP2P1.o \
../boundary_value.o \
../build_fem_space.o \
../build_matrix.o \
../build_matrix_struct.o \
../build_mesh.o \
../config.o \
../controler.o \
../functions.o \
../main.o \
../moving.o \
../output_tecplot.o \
../preconditioner.o \
../solver.o \
../step_forward_Euler.o \
../step_forward_linearized_Euler.o 

CPP_SRCS += \
../ISOP2P1.cpp \
../boundary_value.cpp \
../build_fem_space.cpp \
../build_matrix.cpp \
../build_matrix_struct.cpp \
../build_mesh.cpp \
../config.cpp \
../controler.cpp \
../functions.cpp \
../main.cpp \
../moving.cpp \
../output_tecplot.cpp \
../preconditioner.cpp \
../solver.cpp \
../step_forward_Euler.cpp \
../step_forward_linearized_Euler.cpp 

OBJS += \
./ISOP2P1.o \
./boundary_value.o \
./build_fem_space.o \
./build_matrix.o \
./build_matrix_struct.o \
./build_mesh.o \
./config.o \
./controler.o \
./functions.o \
./main.o \
./moving.o \
./output_tecplot.o \
./preconditioner.o \
./solver.o \
./step_forward_Euler.o \
./step_forward_linearized_Euler.o 

CPP_DEPS += \
./ISOP2P1.d \
./boundary_value.d \
./build_fem_space.d \
./build_matrix.d \
./build_matrix_struct.d \
./build_mesh.d \
./config.d \
./controler.d \
./functions.d \
./main.d \
./moving.d \
./output_tecplot.d \
./preconditioner.d \
./solver.d \
./step_forward_Euler.d \
./step_forward_linearized_Euler.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


