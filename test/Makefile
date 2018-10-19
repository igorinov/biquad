CFLAGS += -I.. -L..

test: test.c ../biquad.o
	$(CC) $(CFLAGS) -o $@ test.c ../biquad.o -lm

clean:
	rm -f test

