CXX ?= g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Iinclude -IImgLib/include -IImgLib/third_party

BUILD_DIR = build
IMGLIB_A = ImgLib/build/libimglib.a
LIB_TARGET = $(BUILD_DIR)/libinpaint.a
CLI_TARGET = $(BUILD_DIR)/inpaint
TEST_TARGET = $(BUILD_DIR)/inpaint_tests

LIB_SOURCES = \
	src/operators.cpp \
	src/deconv.cpp \
	src/poisson.cpp \
	src/inpaint.cpp

LIB_OBJECTS = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(LIB_SOURCES))
CLI_OBJECT = $(BUILD_DIR)/src/cli.o
TEST_OBJECT = $(BUILD_DIR)/tests/test_inpaint.o

all: $(LIB_TARGET) $(CLI_TARGET)

$(IMGLIB_A):
	$(MAKE) -C ImgLib build/libimglib.a

$(LIB_TARGET): $(LIB_OBJECTS)
	ar rcs $@ $^

$(CLI_TARGET): $(CLI_OBJECT) $(LIB_TARGET) $(IMGLIB_A)
	$(CXX) $(CXXFLAGS) -o $@ $(CLI_OBJECT) $(LIB_TARGET) $(IMGLIB_A)

$(TEST_TARGET): $(TEST_OBJECT) $(LIB_TARGET) $(IMGLIB_A)
	$(CXX) $(CXXFLAGS) -o $@ $(TEST_OBJECT) $(LIB_TARGET) $(IMGLIB_A)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

test: $(TEST_TARGET)
	./$(TEST_TARGET)

check: test

clean:
	rm -rf $(BUILD_DIR)
	$(MAKE) -C ImgLib clean

.PHONY: all test check clean
