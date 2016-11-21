CFLAGS = -g -pg -Wall -W -Wextra -O3 --fast-math
all: prog.exe
prog.exe: ./objs/main.o ./objs/heat_pchm.o
	g++ $(CFLAGS) ./objs/main.o ./objs/heat_pchm.o -o prog.exe -lm
./objs/heat_pchm.o: ./heat_model11/heat_pchm.cpp ./heat_model11/func_header.h
	g++ $(CFLAGS) -c ./heat_model11/heat_pchm.cpp -o ./objs/heat_pchm.o
./objs/main.o: ./heat_model11/main.cpp ./heat_model11/func_header.h
	g++ $(CFLAGS) -c ./heat_model11/main.cpp -o ./objs/main.o
clean:
	rm -f ./objs/*.o prog.exe
