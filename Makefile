main: sidh32.c
	clang -c -emit-llvm -o sidh32.bc sidh32.c -Wall -Wextra -Wpedantic -Werror
	saw sidh32.saw --sim-verbose=3

build: sidh32.c
	clang -g -O0 -o sidh32.out sidh32.c -Wall -Wextra -Wpedantic -Werror
	./sidh32.out

clean:
	rm *.bc
	rm *.out
