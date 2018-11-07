CC = g++
EXE = calc
OBJ = flat.o

$(EXE):$(OBJ)
	$(CC) -o $(EXE) -L/usr/local/lib $(OBJ) -lgsl -lgslcblas -lm

%.o:%.c
	$(CC) -wall -I/usr/lcoal/include -c $^ -o $@ 

clean:
	rm -rf $(OBJ) $(EXE)
