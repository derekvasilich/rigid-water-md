CC=gcc
NAME=water
BIN=${NAME}.exe

all:	list.o
	${CC} -o ${BIN} list.o -lGL -lGLU -lglut ${NAME}.c 