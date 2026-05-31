CC = gcc
CFLAGS = -O3 -march=skylake -Wall -Wextra
LDFLAGS = -lpthread -lm

TARGETS = krakken_scalar krakken_multi verify_krakken

.PHONY: all clean run

all: $(TARGETS)

krakken_scalar: krakken.c krakken.h
	$(CC) $(CFLAGS) -DKRAKKEN_MAIN $< -o $@ $(LDFLAGS)

krakken_multi: krakken_multi.c
	$(CC) $(CFLAGS) -DKRAKKEN_MAIN $< -o $@ $(LDFLAGS)

verify_krakken: verify_krakken.c krakken.c krakken_multi.c krakken.h
	$(CC) $(CFLAGS) verify_krakken.c krakken.c krakken_multi.c -o $@ $(LDFLAGS)

run: all
	@echo "=== Running Krakken Verification Tests ==="
	./verify_krakken
	@echo ""
	@echo "=== Running Krakken Scalar Benchmark ==="
	./krakken_scalar
	@echo ""
	@echo "=== Running Krakken AVX2/Parallel Benchmark ==="
	./krakken_multi

clean:
	rm -f $(TARGETS)
