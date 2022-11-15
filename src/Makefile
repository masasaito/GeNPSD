
INC_COMMON = ./Common.mk
include $(INC_COMMON)

## Macros
TARGETS  = genpsd_generate genpsd_interconv
SHLS     = $(STARGETS)
BUPIDX   = $(shell date "+%Y%m%d_%H%M")
BUPFILE  = genpsd-src.$(BUPIDX).tar.gz

## Files
SRCS = main_genpsd_generate.f90 main_genpsd_interconv.f90 genpsd_module.f90
OBJS = $(SRCS:.f90=.o)

## Commands
all : $(TARGETS) $(STARGETS)

install : $(TARGETS) $(STARGETS)
	@ mkdir -p $(BINDIR)
	cp $(TARGETS) $(STARGETS) $(BINDIR)
clean  :
	rm -f *.o *.mod $(TARGETS) *~
echo   :
	@ echo $(SRCS) $(STARGETS)
print  :
	@ cat $(SRCS) $(STARGETS)
backup : $(BUPFILE)
$(BUPFILE) :
	@ rm -f $(BUPFILE)
	tar cvzf $(BUPFILE) $(SRCS) $(STARGETS) Makefile *.mk

## Suffix rules
.f90.o :
	$(FC) -c $(FCFLAGS) $(LIB_FLAGS) $(LAPACK_FLAGS) $<

## Targets
genpsd_generate      : main_genpsd_generate.o genpsd_module.o $(LIB_OBJS) 
	$(FC) -o $@ $(FCFLAGS) $^
genpsd_interconv     : main_genpsd_interconv.o genpsd_module.o $(LIB_OBJS) 
	$(FC) -o $@ $(FCFLAGS) $^

## Dependency
main_genpsd_generate.o  : $(LIB_OBJS) genpsd_module.o
main_genpsd_interconv.o : $(LIB_OBJS) genpsd_module.o
genpsd_module.o         : $(LIB_OBJS)

