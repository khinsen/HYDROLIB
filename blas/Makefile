.f.o:
	$(FC) $(FCFLAGS1) $*.f

OBJS = \
  daxpy.o \
  dgemm.o \
  dnrm2.o \
  dspr.o \
  dsyr.o \
  dsyrk.o \
  dtrmv.o \
  lsame.o  \
  dcopy.o \
  dgemv.o \
  drot.o \
  dswap.o \
  dsyr2.o \
  dtpsv.o \
  dtrsm.o \
  xerbla.o  \
  ddot.o \
  dger.o \
  dscal.o \
  dsymv.o \
  dsyr2k.o \
  dtrmm.o \
  idamax.o \
  saxpy.o \
  sgemm.o \
  snrm2.o \
  sspr.o \
  ssyr.o \
  ssyrk.o \
  strmv.o \
  scopy.o \
  sgemv.o \
  srot.o \
  sswap.o \
  ssyr2.o \
  stpsv.o \
  strsm.o \
  sdot.o \
  sger.o \
  sscal.o \
  ssymv.o \
  ssyr2k.o \
  strmm.o \
  isamax.o

libblas.a: ${OBJS}
	ar r libblas.a ${OBJS}
	$(RANLIB) libblas.a
