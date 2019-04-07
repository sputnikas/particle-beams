CC = g++
CFLAGS = -c -Wall -O0
LDFLAGS = 
DEFINES = -DHAS_TEST #-DHAS_TEST_CONSTRUCTOR
SRC = main.cpp src/vec2.cpp src/vec3.cpp
HDR = $(SRC:.cpp=.h)
OBJ = $(SRC:.cpp=.o)
OBJDIR = obj/
EXE = bin/test.exe
ASM = $(SRC:.cpp=.s)
ASMDIR = asm/

all: dir $(EXE)

dir:
	mkdir -p obj/
	mkdir -p obj/src/

$(EXE): $(addprefix $(OBJDIR), $(OBJ))
	$(CC) $(addprefix $(OBJDIR), $(OBJ)) -o $@ $(LDFLAGS)
	
$(OBJDIR)%.o: %.cpp %.h
	$(CC) $(CFLAGS) $(DEFINES) $< -o $@ 

clean:
	rm -f $(addprefix $(OBJDIR), $(OBJ)) $(EXE)
	

asm : asmdir $(addprefix $(ASMDIR), $(ASM))
	
asmdir : 
	mkdir -p asm/
	mkdir -p asm/src/

%.s: %.cpp
	$(CC) -S -masm=intel $(DEFINES) $< -o $(ASMDIR)$@

cleanasm:
	rm -f $(addprefix $(ASMDIR), $(ASM))