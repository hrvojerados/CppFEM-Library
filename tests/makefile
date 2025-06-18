FLAGS = -O3 -I../thirdParty/include -std=c++17
TARGETS = test1 test2 test3 meshTest

all: $(TARGETS)

test1:
	g++ $(FLAGS) test1.cpp -o test1

test2:
	g++ $(FLAGS) test2.cpp -o test2

test3:
	g++ $(FLAGS) test3.cpp -o test3

meshTest:
	g++ $(FLAGS) meshTest.cpp -o meshTest

clean:
	rm test1
	rm test2
	rm test3
	rm meshTest
