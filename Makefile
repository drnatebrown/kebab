CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native -flto \
           -ffast-math \
           -funroll-loops \
           -fomit-frame-pointer \
           -DNDEBUG
LDFLAGS = -flto -Wl,-O3
CXXFLAGS_DEBUG = -std=c++17 -Wall -Wextra -O0 -g -fsanitize=address,undefined
LDFLAGS_DEBUG = -fsanitize=address,undefined
INCLUDES = -I./include

SRC_DIR = src
OBJ_DIR = obj

SRCS = src/kebab.cpp \
       src/kebab/kebab_index.cpp \
       src/kebab/nt_hash.cpp \
       src/external/hll/hll.cpp
OBJS = obj/kebab.o \
       obj/kebab/kebab_index.o \
       obj/kebab/nt_hash.o \
       obj/external/hll/hll.o

# Add header dependencies
DEPS = $(OBJS:.o=.d)

TARGET = kebab

.PHONY: all clean debug

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

-include $(DEPS)

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

debug: clean
debug: CXXFLAGS = $(CXXFLAGS_DEBUG)
debug: LDFLAGS = $(LDFLAGS_DEBUG)
debug: all