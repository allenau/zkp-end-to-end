OS := $(shell uname -s)
CLANG_VERSION ?= $(shell clang --version | grep version | sed 's/[^0-9]*\([0-9]\{1,\}\.[0-9]\{1,\}\.[0-9]\{1,\}\).*/\1/')
ROOT_DIR ?= ../..
SRC_DIR ?= .
PROJ_NAME ?= libcryptoplus
TEST_NAME ?= cryptoplustest
PROJ_OUT ?= $(ROOT_DIR)/build/cpp
TEST_OUT ?= $(SRC_DIR)/build
COVERAGE_OUT ?= $(TEST_OUT)/coverage

GTEST_DIR:= $(ROOT_DIR)/deps/googletest/googletest
CRYPTOPP_DIR ?= $(ROOT_DIR)/build/cpp/vendor/cryptopp
DJINNI_GEN_DIR ?= $(ROOT_DIR)/generated-src/cpp
CPP_DIR := $(ROOT_DIR)/src/cpp

CPP_SRC_FILES := $(wildcard \
	${CPP_DIR}/*.cpp \
	${CPP_DIR}/**/*.cpp \
	${CPP_DIR}/**/**/*.cpp \
)
CPP_OBJS := $(patsubst ${CPP_DIR}/%,%,${CPP_SRC_FILES})
CPP_INCLUDES := \
	$(CRYPTOPP_DIR) \
	$(DJINNI_GEN_DIR) \
	$(CPP_DIR)
CPP_TARGET := $(wildcard \
	${CPP_DIR}/*.o \
	${CPP_DIR}/**/*.o \
	${CPP_DIR}/**/**/*.o \
)

TEST_FILES := $(wildcard ${SRC_DIR}/*.cpp)
TEST_LIST := $(patsubst ${SRC_DIR}/%.cpp,%,${TEST_FILES})

FLAGS := -std=c++11 -isystem ${GTEST_DIR}/include -pthread
COVERAGE_FLAG:=
ifeq ($(OS),Darwin)
	COVERAGE_FLAG = -L/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/$(CLANG_VERSION)/lib/darwin \
	-lclang_rt.profile_osx
else
	COVERAGE_FLAG = -lgcov
endif

all: compile_program clean compile_test

clean:
	@echo "\033[33mClean unit test project...\033[0m"
	rm -rf $(TEST_OUT)

compile_gtest:
	@echo "\033[33mBuild google test...\033[0m"
	@cd $(GTEST_DIR) && mkdir -p build
	g++ $(FLAGS) -I$(GTEST_DIR) -c $(GTEST_DIR)/src/gtest-all.cc \
		-o $(GTEST_DIR)/build/gtest-all.o
	g++ $(FLAGS) -I$(GTEST_DIR) -c $(GTEST_DIR)/src/gtest_main.cc \
		-o $(GTEST_DIR)/build/gtest_main.o
	ar rv $(GTEST_DIR)/build/gtest_main.a \
		$(GTEST_DIR)/build/gtest-all.o $(GTEST_DIR)/build/gtest_main.o

compile_program:
	@echo "\033[33mBuild core project...\033[0m"
	@cd $(CPP_DIR) && make clean library

compile_test:
	@echo "\033[33mBuild test...\033[0m"
	mkdir -p $(TEST_OUT)

	$(foreach f, $(TEST_LIST), \
		g++ $(FLAGS) \
			$(foreach d, $(CPP_INCLUDES), -I$d) \
			-c $(SRC_DIR)/$f.cpp -o $(TEST_OUT)/$f.o; )

	@echo "\033[33mRun test...\033[0m"
	g++ $(FLAGS) -lpthread $(CPP_TARGET) \
		$(foreach f, $(TEST_LIST), $(TEST_OUT)/$f.o) \
		$(GTEST_DIR)/build/gtest_main.a -o $(TEST_OUT)/$(TEST_NAME) \
			$(COVERAGE_FLAG) \
			-L$(CRYPTOPP_DIR) \
			-lcryptopp;

	@echo "\033[32mOutput:\033[0m\n$(TEST_OUT)"
	@echo "\033[32mRunning:\033[0m"

	$(TEST_OUT)/$(TEST_NAME) --gtest_output="xml:$(TEST_OUT)/$(TEST_NAME).xml";

test_coverage:
	@echo "\033[33mRun coverage report...\033[0m"
	cd $(CPP_DIR) && gcov $(foreach f, $(CPP_OBJS), $f) > /dev/null

	@mkdir -p $(COVERAGE_OUT)
	@lcov -q --no-external -c -b $(CPP_DIR)/lib -d $(CPP_DIR)/lib --ignore-errors graph -o $(COVERAGE_OUT)/coverage.info
	@lcov -r $(COVERAGE_OUT)/coverage.info '$(CPP_DIR)/lib/*/*.hpp' '$(CPP_DIR)/lib/*.hpp' '*.hpp' -o $(COVERAGE_OUT)/coverage.info
	@genhtml $(COVERAGE_OUT)/coverage.info -o $(COVERAGE_OUT)