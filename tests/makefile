FLAGS = -O3 -I../thirdParty/include -std=c++17
TARGETS = test1 test2 test3 test4 meshTest

all: $(TARGETS)

test1:
	g++ $(FLAGS) test1.cpp -o test1

test2:
	g++ $(FLAGS) test2.cpp -o test2

test3:
	g++ $(FLAGS) test3.cpp -o test3

test4:
	g++ $(FLAGS) test4.cpp -o test4

meshTest:
	g++ $(FLAGS) meshTest.cpp -o meshTest

clean:
	rm test1
	rm test2
	rm test3
	rm test4
	rm meshTest
