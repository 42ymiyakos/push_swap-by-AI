NAME    = push_swap
CC      ?= cc
CFLAGS  ?= -O2 -std=gnu11 -Wall -Wextra
SRCS    = main.c util.c rep.c construct.c beam.c exact.c peep.c ig.c
OBJS    = $(SRCS:.c=.o)

all: $(NAME)
$(NAME): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS)
%.o: %.c ps.h
	$(CC) $(CFLAGS) -c $< -o $@
asan:
	$(CC) -O1 -g -std=gnu11 -fsanitize=address,undefined -o $(NAME)_asan $(SRCS)
clean:
	rm -f $(OBJS)
fclean: clean
	rm -f $(NAME) $(NAME)_asan
re: fclean all
.PHONY: all asan clean fclean re
