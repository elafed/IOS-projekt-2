proj2: proj2.o
	gcc	-std=gnu99 -Wall -Wextra -Werror -pedantic -pthread -o proj2 proj2.o

proj2.o: proj2.c
	gcc	-std=gnu99 -Wall -Wextra -Werror -pedantic -pthread -c proj2.c

.PHONY: clean
clean:
	rm proj2 *.o

