# Makefile for Hydrodynamic Interaction library

LIBOBJ = gcalc0.o gcalcpbc0.o gcalc1.o gcalcpbc1.o\
         gcalc2.o gcalcpbc2.o gcalc3.o gcalcpbc3.o\
         lubrication.o lubricationpbc.o

lib: $(LIBOBJ)

gcalc0.o: gcalc.f filenames.h
	./fcomp gcalc$(DPSUFFIX).f $@ -D_LM_=0 $(DPOPT)
gcalc1.o: gcalc.f filenames.h
	./fcomp gcalc$(DPSUFFIX).f $@ -D_LM_=1 $(DPOPT)
gcalc2.o: gcalc.f filenames.h
	./fcomp gcalc$(DPSUFFIX).f $@ -D_LM_=2 $(DPOPT)
gcalc3.o: gcalc.f filenames.h
	./fcomp gcalc$(DPSUFFIX).f $@ -D_LM_=3 $(DPOPT)

gcalcpbc0.o: gcalc.f filenames.h
	./fcomp gcalc$(DPSUFFIX).f $@ -D_LM_=0 -DPERIODIC $(DPOPT)
gcalcpbc1.o: gcalc.f filenames.h
	./fcomp gcalc$(DPSUFFIX).f $@ -D_LM_=1 -DPERIODIC $(DPOPT)
gcalcpbc2.o: gcalc.f filenames.h
	./fcomp gcalc$(DPSUFFIX).f $@ -D_LM_=2 -DPERIODIC $(DPOPT)
gcalcpbc3.o: gcalc.f filenames.h
	./fcomp gcalc$(DPSUFFIX).f $@ -D_LM_=3 -DPERIODIC $(DPOPT)

lubrication.o: lubrication.f filenames.h
	./fcomp lubrication$(DPSUFFIX).f $@ $(DPOPT)
lubricationpbc.o: lubrication.f filenames.h
	./fcomp lubrication$(DPSUFFIX).f $@ -DPERIODIC $(DPOPT)
