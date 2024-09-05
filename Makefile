# Build environment can be configured the following
# environment variables:
#   CC : Specify the C compiler to use
#   CFLAGS : Specify compiler options to use
#   LDFLAGS : Specify linker options to use
#   CPPFLAGS : Specify c-preprocessor options to use
CC = gcc
CFLAGS = -Wall -O3 -g -flto

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

all: clean static

all_fma: CFLAGS+=-mfma
all_fma: CFLAGS+=-D__USE_FMA_INTRIN
all_fma: clean static

all_avx: CFLAGS+=-mfma
all_avx: CFLAGS+=-mavx
all_avx: CFLAGS+=-D__USE_AVX_INTRIN
all_avx: CFLAGS+=-D__USE_FMA_INTRIN
all_avx: clean static

static: $(LIB_A)

shared dynamic: $(LIB_SO)

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

clean:
	@echo "Cleaning build objects & library..."
	@$(RM) $(LIB_OBJS) $(LIB_LOBJS) $(LIB_A) $(LIB_SO) $(LIB_SO_MAJOR) $(LIB_SO_BASE)
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
