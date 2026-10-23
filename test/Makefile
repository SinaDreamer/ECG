CC      ?= cc
CFLAGS  ?= -std=c99 -O2 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion
DSP     := ../dsp
SRCS    := test_ecg.c $(wildcard $(DSP)/*.c)

test_ecg: $(SRCS) $(wildcard $(DSP)/*.h)
	$(CC) $(CFLAGS) -I$(DSP) -o $@ $(SRCS) -lm

run: test_ecg
	./test_ecg

clean:
	rm -f test_ecg

.PHONY: run clean
