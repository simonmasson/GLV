CC=gcc
CFLAGS= -lgmp -lm
DEPS = entete.h
OBJ =  main.o points.o add_points.o multiple_point.o double_and_add.o scalar_decomposition.o simultaneous_multiple_multiplication.o signed_binary.o jacobian_to_weierstrass.o
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	gcc $^ $(CFLAGS)
