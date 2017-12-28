################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
main.cpp 

CC_SRCS += \
io/batch.cc \
fea/fb.cc \
fea/fea.cc \
fea/post.cc \
fea/fea_impl.cc \
fea/post_impl.cc \
fea/fea_trap.cc \
fea/fea_delta.cc \
io/in.cc \
nr/nr.cc \
io/opts.cc \
io/out.cc \
io/pfile.cc \
vad/vad.cc

CPP_DEPS += \
main.d 

CC_DEPS += \
io/batch.d \
fea/fb.d \
fea/fea.d \
fea/post.d \
fea/fea_impl.d \
fea/post_impl.d \
fea/fea_trap.d \
fea/fea_delta.d \
io/in.d \
nr/nr.d \
io/opts.d \
io/out.d \
io/pfile.d \
vad/vad.d

OBJS += \
io/batch.o \
fea/fb.o \
fea/fea.o \
fea/post.o \
fea/fea_impl.o \
fea/post_impl.o \
fea/fea_trap.o \
fea/fea_delta.o \
io/in.o \
main.o \
nr/nr.o \
io/opts.o \
io/out.o \
io/pfile.o \
vad/vad.o


# Each subdirectory must supply rules for building sources it contributes
%.o: %.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: %.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


