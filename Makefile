all:
	rm -f App
	g++ -std=c++11 basic_primitives.cpp he_lbp.cpp main.cpp ./../fhe.a -o App -lntl -lgmp -lm
	./App L=30 s=8 nt=2
clean:
	rm -f App
