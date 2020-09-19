CC = gcc
FLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LIBS = -lm
OBJFILES = cluster.o spmat.o Algorithm3.o error_codes.o modularity_maximization.o Algorithm2.o  one_norm.o power_iter.o 
TARGET = cluster

all: $(TARGET)
$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES)${LIBS}

.PHONY: clean
clean:
	rm -f $(OBJFILES) $(TARGET) *~

