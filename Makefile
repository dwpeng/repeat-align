CC = gcc
CFLAGS = -Wall -Wextra -pedantic -O0 -g

all: repeat-align

repeat-align: repeat-align.c
	$(CC) $(CFLAGS) -o repeat-align repeat-align.c

clean:
	rm -f repeat-align
