all:
	rm -f App
	g++ -std=c++11 basic_primitives.cpp main.cpp ./../fhe.a -o App -lntl -lgmp -lm
	./App L=40 s=1024
clean:
	rm -f App
