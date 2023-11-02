CXX = g++
CXXFLAGS = -std=c++11 -Iinclude -IImgLib/include
SRC_DIR = src
BUILD_DIR = build
EXECUTABLE = imglib

# List of source files including those from the submodule
SOURCES := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard ImgLib/src/*.cpp)

# List of object files (one for each source file)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))

all: $(BUILD_DIR)/$(EXECUTABLE)

$(BUILD_DIR)/$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Pattern rule to build object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
