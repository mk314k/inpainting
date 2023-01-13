NUMBER := 9
OBJS=\
$(BUILD_DIR)/basicImageManipulation.o \
$(BUILD_DIR)/a9.o \
$(BUILD_DIR)/morphing.o \
$(BUILD_DIR)/filtering.o \
$(BUILD_DIR)/npr.o \
$(BUILD_DIR)/Image.o \
$(BUILD_DIR)/lodepng.o 

include Makefile.include

prepare:
	mkdir -p asst
	cp *.cpp asst
	cp *.h asst
	cp *.py asst
	cp Makefile asst
	cp Makefile.include asst
	cp -R Input asst
	cp -R Eigen asst
	cp -R Output asst
	cp -R write-up asst
	zip -r a$(NUMBER)_submission.zip asst
	rm -rf asst
.PHONY: prepare