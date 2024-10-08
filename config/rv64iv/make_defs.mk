#
#
#  BLIS
#  An object-based framework for developing high-performance BLAS-like
#  libraries.
#
#  Copyright (C) 2014, The University of Texas at Austin
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are
#  met:
#   - Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   - Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   - Neither the name(s) of the copyright holder(s) nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#


# Declare the name of the current configuration and add it to the
# running list of configurations included by common.mk.
THIS_CONFIG    := rv64iv
#CONFIGS_INCL   += $(THIS_CONFIG)

#
# --- Determine the C compiler and related flags ---
#

# NOTE: The build system will append these variables with various
# general-purpose/configuration-agnostic flags in common.mk. You
# may specify additional flags here as needed.
CPPROCFLAGS    := -DRISCV_SIZE=64

#RISCV_ARCH     := $(shell $(CC) -DFORCE_RISCV_VECTOR -E build/detect/riscv/bli_riscv_detect_arch.h | grep '^[^\#]')
#RISCV_ABI      := $(shell $(CC) -DFORCE_RISCV_VECTOR -E build/detect/riscv/bli_riscv_detect_abi.h | grep '^[^\#]')
RISCV_ARCH	:= rv64gc_zfh_xtheadvector
RISCV_ABI	:= lp64d

ifeq (,$(findstring 64,$(RISCV_ARCH)))
$(error The RISC-V compiler architecture $(RISCV_ARCH) is not compatible with $(THIS_CONFIG))
else ifeq (,$(findstring 64,$(RISCV_ABI)))
$(error The RISC-V compiler ABI $(RISCV_ABI) is not compatible with $(THIS_CONFIG))
endif

CMISCFLAGS     := -march=$(RISCV_ARCH) -mabi=$(RISCV_ABI)
CPICFLAGS      := -fPIC
CWARNFLAGS     := -Wall -Wno-unused-function -Wfatal-errors

# In case the A extension is not available
LDFLAGS        += -latomic

ifneq ($(DEBUG_TYPE),off)
CDBGFLAGS      := -g
endif

ifeq ($(DEBUG_TYPE),noopt)
COPTFLAGS      := -O0
else
COPTFLAGS      := -O2 -ftree-vectorize
endif

# Flags specific to optimized kernels.
CKOPTFLAGS     := $(COPTFLAGS) -O3
ifeq ($(CC_VENDOR),gcc)
CKVECFLAGS     :=
else
ifeq ($(CC_VENDOR),clang)
CKVECFLAGS     :=
else
$(error gcc or clang is required for this configuration.)
endif
endif

# Flags specific to reference kernels.
CROPTFLAGS     := $(CKOPTFLAGS)
ifeq ($(CC_VENDOR),gcc)
# Lower compiler optimization. cinvscalv fails at -O1
CRVECFLAGS     := $(CKVECFLAGS) -O0
else
ifeq ($(CC_VENDOR),clang)
CRVECFLAGS     := $(CKVECFLAGS) -funsafe-math-optimizations -ffp-contract=fast
else
CRVECFLAGS     := $(CKVECFLAGS)
endif
endif

# Store all of the variables here to new variables containing the
# configuration name.
$(eval $(call store-make-defs,$(THIS_CONFIG)))
