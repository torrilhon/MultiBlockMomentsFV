# ==========================================================================
#  MultiBlockMomentsFV  --  build / run / check / clean
#
#  Quick start:
#     make                         build build/moments
#     make run CASE=cases/squareSmall.txt
#     make check                   regression run vs reference/ (see note below)
#     make clean                   remove build artifacts
#
#  (c) Manuel Torrilhon, RWTH Aachen 2026
# ==========================================================================

CC      := gcc
CFLAGS  := -O3 -Wall
LDFLAGS := -lm

SRCDIR  := src
BUILD   := build
RUNS    := runs
REF     := reference

# all translation units
SRC     := $(wildcard $(SRCDIR)/*.c)
OBJ     := $(patsubst $(SRCDIR)/%.c,$(BUILD)/%.o,$(SRC))
BIN     := $(BUILD)/moments

# default case for `make run`
CASE    ?= cases/squareSmall.txt
# default case for `make check` (square4 = full standard reference run)
CHECKCASE ?= cases/square4.txt
CHECKNAME  = $(basename $(notdir $(CHECKCASE)))
# tolerances for the numeric regression comparison
ABSTOL  ?= 1e-10
RELTOL  ?= 1e-8
# strip directory + .txt extension -> name of the run/output subfolder
CASENAME = $(basename $(notdir $(CASE)))

.PHONY: all run check clean help

# --------------------------------------------------------------------------
all: $(BIN)

$(BIN): $(OBJ) | $(BUILD)
	$(CC) $(OBJ) $(LDFLAGS) -o $@
	@echo "  -> built $(BIN)"

$(BUILD)/%.o: $(SRCDIR)/%.c | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD):
	@mkdir -p $(BUILD)

# --------------------------------------------------------------------------
# run: each case writes into its own runs/<case>/ directory via the -out flag.
# Pass extra solver flags through ARGS, e.g.  make run CASE=cases/square4.txt ARGS="-Kn 0.01"
run: $(BIN)
	@mkdir -p $(RUNS)/$(CASENAME)
	./$(BIN) -bfile $(CASE) -out $(RUNS)/$(CASENAME) $(ARGS)

# --------------------------------------------------------------------------
# check: run CHECKCASE and compare its Tecplot output against reference/
# numerically (tools/compare_solution.py). The comparison ignores physically
# meaningless differences (signed zeros, whitespace, Tecplot STRANDID) and
# fails on any value that exceeds the tolerance. This is a hard pass/fail gate.
#
#   make check                          # default: square4
#   make check CHECKCASE=cases/squareSmall.txt   # quick smoke check
#   make check ABSTOL=1e-8 RELTOL=1e-6           # loosen tolerances
check: $(BIN)
	@echo "== regression check on $(CHECKCASE) =="
	@if [ ! -f $(REF)/$(CHECKNAME).ref.dat ]; then \
	   echo "SKIP: no reference file $(REF)/$(CHECKNAME).ref.dat"; exit 0; \
	 fi
	@mkdir -p $(RUNS)/$(CHECKNAME)
	./$(BIN) -bfile $(CHECKCASE) -out $(RUNS)/$(CHECKNAME)
	@echo "-- compare vs $(REF)/$(CHECKNAME).ref.dat --"
	python3 tools/compare_solution.py \
	   $(RUNS)/$(CHECKNAME)/*_tAll.dat $(REF)/$(CHECKNAME).ref.dat $(ABSTOL) $(RELTOL)

# --------------------------------------------------------------------------
clean:
	rm -rf $(BUILD)
	@echo "  -> removed $(BUILD)/  (run results in $(RUNS)/ are kept)"

help:
	@echo "targets: all (default), run, check, clean"
	@echo "vars:    CASE=cases/<file>.txt        ARGS=\"-Kn 0.01 -CFL 0.9\""
	@echo "         CHECKCASE=cases/<file>.txt   ABSTOL=1e-10 RELTOL=1e-8"
