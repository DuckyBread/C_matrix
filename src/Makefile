CC = gcc
CFLAGS =-Wall -Wextra -Werror -g -std=c11
CFLAGS_GCOV =-fprofile-arcs --coverage

OS=$(shell uname)
ifeq ($(OS), Linux)
  LIBS=-lcheck -lsubunit -lpthread -lrt -lm -g
  CHECK_LEAKS=valgrind --leak-check=full -s --log-file=log.txt
  OPEN_GCOV=firefox coverage/src/index.html
else
  LIBS=-lcheck 
  CHECK_LEAKS=no leaks -atExit -- 
  OPEN_GCOV=open coverage/index.html
endif

all: my_matrix.a test clean

test: my_matrix.o my_matrix.a
	checkmk clean_mode=1 my_matrix_test.check > my_matrix_test.c
	gcc $(CFLAGS_GCOV) my_matrix_test.c my_matrix.c $(LIBS)
	./a.out

my_matrix.a:
	$(CC) $(CFLAGS) -c my_matrix.c -o my_matrix.o
	ar rsc my_matrix.a my_matrix.o
	
style:
	cp ../materials/linters/.clang-format .
	clang-format -i *.c
	clang-format -n *.c
	rm .clang-format

gcov_report:
	gcovr --html-details coverage.html
	open coverage.html

leaks: test
	$(CHECK_LEAKS) ./a.out

clean:
	rm -rf ./lcov_report ./*.dSYM
	rm -f *.o *.a *.gcno *.gcda *.html *.css *.out *.gch *.cfg coverage *.out *.txt my_matrix_test.c 
