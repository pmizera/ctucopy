################################################################################
# Automatically-generated file. Do not edit!
################################################################################

RM := rm

# All of the sources participating in the build are defined here
-include sources.mk
-include vdet/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CU_DEPS)),)
-include $(CU_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(TCC_DEPS)),)
-include $(TCC_DEPS)
endif
endif

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: ctucopy4

# Tool invocations
ctucopy4: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++  -o "bin/ctucopy4" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(C_DEPS)$(C++_DEPS)$(CXX_DEPS)$(CPP_DEPS)$(CU_DEPS)$(CC_DEPS)$(EXECUTABLES)$(OBJS)$(C_UPPER_DEPS)$(TCC_DEPS) bin/ctucopy4
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
