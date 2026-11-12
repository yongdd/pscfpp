include ../../../config.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/pscf/sources.mk
include $(SRC_DIR)/pscf/patterns.mk

TEST=$(BLD_DIR)/pscf/tests/regression/Test

all: $(TEST)

ifdef PSCF_CUDA
TEST_OBJ=$(TEST).ou
TEST_DEP=$(TEST).du

$(TEST): $(TEST).ou $(PSCF_LIBS)
	@cp $(TEST).ou $(TEST).o
	$(NVXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
	@rm $(TEST).o
else
TEST_OBJ=$(TEST).o
TEST_DEP=$(TEST).d

$(TEST): $(TEST).o $(PSCF_LIBS)
	$(CXX) $(LDFLAGS) -o $(TEST) $(TEST).o $(LIBS)
endif

run: $(TEST)
	@$(TEST) $(SRC_DIR)/pscf/tests/regression/

run-quiet: $(TEST_OBJ) $(TEST)
	@$(TEST) $(SRC_DIR)/pscf/tests/regression/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pscf/tests/regression/log" > count
	@cat count

clean-outputs:
	rm -f log count

clean:
	rm -f $(TEST_OBJ)
	rm -f $(TEST_DEP)
	rm -f $(TEST)
	$(MAKE) clean-outputs

.PHONY: all run run-quiet clean-outputs clean

-include $(util_DEPS)
-include $(pscf_DEPS)
-include $(TEST_DEP)
