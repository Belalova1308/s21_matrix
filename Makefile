CC=g++ -std=c++17 -Wall -Werror -Wextra -pedantic -g
FILES= *.cpp
TESTS=tests/*.cpp
CHECKFLAGS=-lgtest -lgtest_main -lpthread
REPORTDIR=gcov_report
GCOV=--coverage

.PHONY: all clean test

all: s21_matrix_oop.a

s21_matrix_oop.a:
	$(CC) -c $(FILES)
	ar rcs s21_matrix_oop.a *.o
	ranlib s21_matrix_oop.a

test: clean
	$(CC) $(GCOV) -c $(FILES) 
	$(CC) -c $(TESTS) 
	$(CC) $(GCOV) -o matrix_test *.o $(CHECKFLAGS)
	./matrix_test

style: 
	cp ../materials/linters/.clang-format .
	clang-format -i *.h $(FILES) $(TESTS)
	rm .clang-format

style_check:
	cp ../materials/linters/.clang-format .
	clang-format -n *.h $(FILES) $(TESTS)
	rm .clang-format

leaks: test
	valgrind -s ./matrix_test

gcov_report: test
	lcov -t "Unit-tests of matrix_oop" -o s21_matrix_oop.info -c -d .
	genhtml -o $(REPORTDIR) s21_matrix_oop.info
	open ./$(REPORTDIR)/index.html

clean:
	rm -rf ./*.o ./*.a ./a.out ./*.gcno ./*.gcda ./$(REPORTDIR) *.info ./*.info report matrix_test matrix_oop
