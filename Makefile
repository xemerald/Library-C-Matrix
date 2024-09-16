# Build environment can be configured the following
# environment variables:
#   CC : Specify the C compiler to use
#   CFLAGS : Specify compiler options to use
#   LDFLAGS : Specify linker options to use
#   CPPFLAGS : Specify c-preprocessor options to use
CC = gcc
CFLAGS = -Wall -O3 -g -flto
CPPFLAGS =
INTRIN_FLAGS =

# Extract version from matrix.h, expected line should include LIBMATRIX_VERSION "#.#.#"
MAJOR_VER = $(shell grep LIBMATRIX_VERSION ./include/matrix.h | grep -Eo '[0-9]+.[0-9]+.[0-9]+' | cut -d . -f 1)
FULL_VER = $(shell grep LIBMATRIX_VERSION ./include/matrix.h | grep -Eo '[0-9]+.[0-9]+.[0-9]+')
COMPAT_VER = $(MAJOR_VER).0.0

# Default settings for install target
PREFIX ?= /usr/local
EXEC_PREFIX ?= $(PREFIX)
LIBDIR ?= $(DESTDIR)$(EXEC_PREFIX)/lib
INCLUDEDIR ?= $(DESTDIR)$(PREFIX)/include/libmatrix
DATAROOTDIR ?= $(DESTDIR)$(PREFIX)/share

LIB_SRCS = ./src/matrix.c

LIB_OBJS = $(LIB_SRCS:.c=.o)
LIB_LOBJS = $(LIB_SRCS:.c=.lo)

LIB_NAME = libmatrix
LIB_A = $(LIB_NAME).a

OS := $(shell uname -s)
AVX_FLAG := $(shell lscpu | grep -io ' avx ' | tr -d '[:space:]')
FMA_FLAG := $(shell lscpu | grep -io ' fma ' | tr -d '[:space:]')

UNIT_TEST = matrix_test

# Build dynamic (.dylib) on macOS/Darwin, otherwise shared (.so)
ifeq ($(OS), Darwin)
	LIB_SO_BASE = $(LIB_NAME).dylib
	LIB_SO_MAJOR = $(LIB_NAME).$(MAJOR_VER).dylib
	LIB_SO = $(LIB_NAME).$(FULL_VER).dylib
	LIB_OPTS = -dynamiclib -compatibility_version $(COMPAT_VER) -current_version $(FULL_VER) -install_name $(LIB_SO)
else
	LIB_SO_BASE = $(LIB_NAME).so
	LIB_SO_MAJOR = $(LIB_NAME).so.$(MAJOR_VER)
	LIB_SO = $(LIB_NAME).so.$(FULL_VER)
	LIB_OPTS = -shared -Wl,--version-script=version.map -Wl,-soname,$(LIB_SO_MAJOR)
endif

# Checking for the CPU flags, if there is the specific flag, turn it on.
ifeq ($(FMA_FLAG), fma)
	INTRIN_FLAGS+=-mfma
	INTRIN_FLAGS+=-D__USE_FMA_INTRIN
endif
#
ifeq ($(AVX_FLAG), avx)
	INTRIN_FLAGS+=-mavx
	INTRIN_FLAGS+=-D__USE_AVX_INTRIN
endif

# Building rules
default: CPPFLAGS+=$(INTRIN_FLAGS)
default: clean static

naive: clean static

static: $(LIB_A)

shared dynamic: $(LIB_SO)

test: clean_test static $(UNIT_TEST)
	@./$(UNIT_TEST)
	@$(RM) ./$(UNIT_TEST)

# Build static library
$(LIB_A): $(LIB_OBJS)
	@echo "Building static library $(LIB_A)..."
	@$(RM) $(LIB_A)
	@$(AR) -crs $(LIB_A) $(LIB_OBJS)

# Build shared/dynamic library
$(LIB_SO): $(LIB_LOBJS)
	@echo "Building shared library $(LIB_SO)..."
	@$(RM) $(LIB_SO) $(LIB_SO_MAJOR) $(LIB_SO_BASE)
	@$(CC) $(CFLAGS) $(LIB_OPTS) -o $(LIB_SO) $(LIB_LOBJS)
	@ln -s $(LIB_SO) $(LIB_SO_BASE)
	@ln -s $(LIB_SO) $(LIB_SO_MAJOR)

$(UNIT_TEST): ./test/munit/munit.c ./test/$(UNIT_TEST).c
	@echo "Compiling $@..."
	@$(CC) $(CFLAGS) -o $@ ./test/munit/munit.c ./test/$(UNIT_TEST).c $(LIB_A)

clean:
	@echo "Cleaning build objects & library..."
	@$(RM) $(LIB_OBJS) $(LIB_LOBJS) $(LIB_A) $(LIB_SO) $(LIB_SO_MAJOR) $(LIB_SO_BASE)
	@echo "All clean."

clean_test:
	@echo "Cleaning unit testing file..."
	@$(RM) $(UNIT_TEST)
	@echo "All clean."

install: shared
	@echo "Installing into $(PREFIX)"
	@mkdir -p $(INCLUDEDIR)
	@cp *.h $(INCLUDEDIR)
	@cp -a $(LIB_SO_BASE) $(LIB_SO_MAJOR) $(LIB_SO_NAME) $(LIB_SO) $(LIBDIR)

.SUFFIXES: .c .o .lo

# Standard object building
.c.o:
	@echo "Compiling $<..."
	@$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Standard object building for shared library using -fPIC
.c.lo:
	@echo "Compiling $<..."
	@$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@

# Print Makefile expanded variables, e.g. % make print-LIB_SO
print-%:
	@echo '$*=$($*)'

FORCE:
