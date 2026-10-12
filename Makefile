CC=gcc
CFLAGS=-Wall -pedantic -O3 -Isrc
ifeq ($(OS), Windows_NT)
	EXE=bin/quadratic_module_checker.exe
else
	EXE=bin/quadratic_module_checker
endif

SRC=$(wildcard src/*.c)
OBJ=$(patsubst src/%.c, obj/%.o, $(SRC))

DIRS=bin obj
$(info $(shell mkdir -p $(DIRS)))


all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@

$(OBJ): obj/%.o : src/%.c
	$(CC) $(CFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm -rf obj $(EXE)