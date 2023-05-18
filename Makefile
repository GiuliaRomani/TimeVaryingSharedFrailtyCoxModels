# CPPFLAGS: optimization of level3 is performed, almost all warnings are activated

# The user has to change this. Fill it with the path of the PACS Exampled directory
PACS_PATH?=/home/giuliaromani/Desktop/Giulia_PACS/pacs-examples/Examples

# PACS libraries are stored here
export PACS_LIB_DIR=$(PACS_PATH)/lib

# Main PACS include directory
export PACS_INC_DIR=$(PACS_PATH)/include

XX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -pedantic -I. -I$(PACS_INC_DIR) -I${mkEigenInc}
LDFLAGS  ?= -L$(PACS_LIB_DIR) -Wl,-rpath=$(PACS_LIB_DIR)
LDLIBS	?= -lpacs
LINK.o := $(LINK.cc) # Use C++ linker.

DEPEND = make.dep

EXEC = main
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

.PHONY = all $(EXEC) $(OBJS) clean distclean $(DEPEND)

all: $(DEPEND) $(EXEC)

$(EXEC): $(OBJS)

$(OBJS): %.o: %.cpp

clean:
	$(RM) $(DEPEND)
	$(RM) *.o

distclean: clean
	$(RM) $(EXEC)
	$(RM) *.dat *.png *.out *.bak *~

$(DEPEND): $(SRCS)
	@$(RM) $(DEPEND)
	@for file in $(SRCS); \
	do \
	  $(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM $${file} >> $(DEPEND); \
	done

-include $(DEPEND)