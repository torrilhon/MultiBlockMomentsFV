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

# short git commit hash + committer date baked in at compile time for output
# provenance (both fall back to "unknown" outside a git checkout)
GIT_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
GIT_DATE   := $(shell git show -s --format=%cI HEAD 2>/dev/null || echo unknown)
CFLAGS += -DGIT_COMMIT=\"$(GIT_COMMIT)\" -DGIT_DATE=\"$(GIT_DATE)\"

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
# check: run CHECKCASE and compare its binary output (.hdr/.dat) against
# reference/ numerically (tools/compare_solution.py). The comparison ignores
# physically meaningless metadata (STRANDID / provenance) and fails on any
# value that exceeds the tolerance. This is a hard pass/fail gate.
#
#   make check                          # default: square4
#   make check CHECKCASE=cases/squareSmall.txt   # quick smoke check
#   make check ABSTOL=1e-8 RELTOL=1e-6           # loosen tolerances
check: $(BIN)
	@echo "== regression check on $(CHECKCASE) =="
	@if [ ! -f $(REF)/$(CHECKNAME).ref.hdr ]; then \
	   echo "SKIP: no reference file $(REF)/$(CHECKNAME).ref.hdr"; exit 0; \
	 fi
	@mkdir -p $(RUNS)/$(CHECKNAME)
	./$(BIN) -bfile $(CHECKCASE) -out $(RUNS)/$(CHECKNAME)
	@echo "-- compare vs $(REF)/$(CHECKNAME).ref.hdr --"
	python3 tools/compare_solution.py \
	   $(RUNS)/$(CHECKNAME)/*_tAll.hdr $(REF)/$(CHECKNAME).ref $(ABSTOL) $(RELTOL)

# --------------------------------------------------------------------------
clean:
	rm -rf $(BUILD)
	@echo "  -> removed $(BUILD)/  (run results in $(RUNS)/ are kept)"

help:
	@echo "targets: all (default), run, check, clean"
	@echo "vars:    CASE=cases/<file>.txt        ARGS=\"-Kn 0.01 -CFL 0.9\""
	@echo "         CHECKCASE=cases/<file>.txt   ABSTOL=1e-10 RELTOL=1e-8"
