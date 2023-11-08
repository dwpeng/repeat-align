CC = gcc
CFLAGS = -Wall -Wextra -pedantic -O0 -g

all: repeat-align

repeat-align: repeat-align.c
	$(CC) $(CFLAGS) -o repeat-align repeat-align.c

align:
	./align.sh

clean:
	rm -f repeat-align
