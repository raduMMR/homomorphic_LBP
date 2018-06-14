all:
	rm -f App
	g++ -std=c++11 basic_primitives.cpp he_lbp.cpp HE_FR_LBP.cpp main.cpp ./../fhe.a -o App -lntl -lgmp -lm
	./App L=30 s=1024
clean:
	rm -f App
