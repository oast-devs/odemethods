CC = gcc -Wall -Wextra

SOURCEDIRS = . 
vpath %.c $(SOURCEDIRS)
CFLAGS += ${addprefix -I,$(SOURCEDIRS)} -march=native -O2 -fPIC
VECTFLAGS = -ftree-vectorize -fopt-info-vec-optimized -fopt-info-vec-missed 

srcs = $(wildcard *.c)
objs = $(srcs:.c=.o)
exec = main


main: $(objs)
	$(CC) $^ -shared -o $(exec) $(CFLAGS)

# vectorized function compilation
ODE_methods.o: ODE_methods.c
	$(CC) $(CFLAGS) $(VECTFLAGS) -c $< -o $@

# a generic method for compiling all source .c files into
# object .o files with automatic creation of file dependencies
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

-include $(deps)

# "clean" target, to remove all binaries ("make clean")
.PHONY: clean
clean:
	rm -f *.o *.d $(ALL_TARGETS)
