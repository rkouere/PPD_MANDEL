# Generated automatically from Makefile.in by configure.

##### User configurable options #####

MAKE          = make --no-print-directory

MPI_CC        = /home/enseign/PPD/mpich/bin/mpicc
MPI_F77       = /home/enseign/PPD/mpich/bin/mpif77
MPI_CLINKER   = $(MPI_CC)
MPI_FLINKER   = $(MPI_F77)
#
DEBUG_LIBS    = 
MPI_INC       = 
MPI_LIBS      = 
MPE_CFLAGS    =  -DMPI_LINUX -DUSE_STDARG -DHAVE_PROTOTYPES 
MPE_FFLAGS    = 
MPE_INCDIR    = /home/enseign/PPD/mpich/include
MPE_LIBDIR    = /home/enseign/PPD/mpich/lib
PROF_LIBS     = -L$(MPE_LIBDIR) -lmpe $(MPI_LIBS)
LOG_LIBS      = -L$(MPE_LIBDIR) -llmpe -lmpe $(MPI_LIBS)
TRACE_LIBS    = -L$(MPE_LIBDIR) -ltmpe $(MPI_LIBS)
ANIM_LIBS     = -L$(MPE_LIBDIR) -lampe -lmpe -lX11 $(MPI_LIBS)
XGRAPH_LIBS   = -L$(MPE_LIBDIR) -lmpe -lX11 $(MPI_LIBS)
F2CMPI_LIBS   = -lfmpich
FLIB_PATH     = -L$(MPE_LIBDIR) $(F2CMPI_LIBS) -llmpe -lmpe \
                $(MPI_LIBS)
#
INCLUDE_DIR   = -I$(MPE_INCDIR) $(MPI_INC)
LDFLAGS       = 
NO_X          = 0
NOF77         = 0
RM            = rm

### End User configurable options ###

VPATH=.:$(srcdir)

srcdir        = /home/enseign/PPD/tmp/mpich-1.2.7p1/mpe/contrib/test

CFLAGS        = $(MPE_CFLAGS) $(INCLUDE_DIR) 
FFLAGS        = $(MPE_FFLAGS) $(INCLUDE_DIR)
CC_LDFLAGS    = $(LDFLAGS) $(MPE_CFLAGS)
FC_LDFLAGS    = $(LDFLAGS) $(MPE_FFLAGS)
EXECS         = cpilog srtest


.c.o:
	$(MPI_CC) $(CFLAGS) -c $<
.f.o:
	$(MPI_F77) $(FFLAGS) -c $<

default: $(EXECS)
	-if [ $(NOF77) = 0 ] ; then \
	     $(MAKE) fpilog ; \
	 fi

cpi_log: cpi.o
	$(MPI_CLINKER) $(CC_LDFLAGS) -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm

cpi_trace: cpi.o
	$(MPI_CLINKER) $(CC_LDFLAGS) -o $@ $? $(TRACE_LIBS) $(DEBUG_LIBS) -lm

cpi_anim: cpi.o
	$(MPI_CLINKER) $(CC_LDFLAGS) -o $@ $? $(ANIM_LIBS) $(DEBUG_LIBS) -lm 

cpilog: cpilog.o
	$(MPI_CLINKER) $(CC_LDFLAGS)  -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm 

fpilog: fpi.o
	$(MPI_FLINKER) $(FC_LDFLAGS) -o $@ $? $(FLIB_PATH) $(DEBUG_LIBS) -lm

srtest: srtest.o
	$(MPI_CLINKER) $(CC_LDFLAGS)  -o $@ $? $(LOG_LIBS) $(DEBUG_LIBS) -lm
cxgraphics: cxgraphics.o
	$(MPI_CLINKER) $(CC_LDFLAGS)  -o $@ $? $(XGRAPH_LIBS) $(DEBUG_LIBS) -lm

fxgraphics: fxgraphics.o
	$(MPI_FLINKER) $(FC_LDFLAGS) -o $@ $? $(XGRAPH_LIBS) $(DEBUG_LIBS) -lm

linktest: linktest_C

linktest_C: 
	@-echo ; \
	echo "** Testing if C application can be linked with logging library" ; \
	$(RM) -f cpi_log ; \
	$(MAKE) cpi_log ; \
	if test -x cpi_log ; then \
	    echo "** C application can be linked with logging library" ; \
	    $(RM) -f cpi_log ; \
	else \
	    echo "** C application CANNOT be linked with logging library" ; \
	fi 
	@-echo ; \
	echo "** Testing if C application can be linked with tracing library" ; \
	$(RM) -f cpi_trace ; \
	$(MAKE) cpi_trace ; \
	if test -x cpi_trace ; then \
	    echo "** C application can be linked with tracing library" ; \
	    $(RM) -f cpi_trace ; \
	else \
	    echo "** C application CANNOT be linked with tracing library" ; \
	fi 
	@-if [ $(NO_X) = 0 ] ; then \
	    echo ; \
	    echo "** Testing if C application can be linked with animation library" ; \
	    $(RM) -f cpi_anim ; \
	    $(MAKE) cpi_anim ; \
	    if test -x cpi_anim ; then \
	        echo "** C application can be linked with animation library" ; \
	        $(RM) -f cpi_anim ; \
	    else \
	        echo "** C application CANNOT be linked with animation library" ; \
	    fi ; \
	    echo ; \
	    echo "** Testing if C application can be linked with graphics library" ; \
	    $(RM) -f cxgraphics ; \
	    $(MAKE) cxgraphics ; \
	    if test -x cxgraphics ; then \
	        echo "** C application can be linked with graphics library" ; \
	        $(RM) -f cxgraphics ; \
	    else \
	        echo "** C application CANNOT be linked with graphics library" ; \
	    fi ; \
	fi  
	@-echo ; \
	echo "** Testing if C application can use both automatic and manual logging together" ; \
	$(RM) -f cpilog ; \
	$(MAKE) cpilog ; \
	if test -x cpilog ; then \
	    echo "** C application can use both automatic and manual logging together" ; \
	    $(RM) -f cpilog ; \
	else \
	    echo "** C application CANNOT use both automatic and manual logging together" ; \
	fi ; \
	echo
	@-$(MAKE) clean

linktest_f77: 
	@-if [ $(NOF77) = 0 ] ; then \
	    echo ; \
	    echo "** Testing if Fortran77 application can be linked with logging library" ; \
	    $(RM) -f fpilog ; \
	    $(MAKE) fpilog ; \
	    if test -x fpilog ; then \
	        echo "** Fortran77 application can be linked with logging library" ; \
	        $(RM) -f fpilog ; \
	    else \
	        echo "** Fortran77 application CANNOT be linked with logging library" ; \
	    fi ; \
	    if [ $(NO_X) = 0 ] ; then \
	        echo ; \
	        echo "** Testing if Fortran77 application can be linked with graphics library" ; \
	        $(RM) -f fxgraphics ; \
	        $(MAKE) fxgraphics ; \
	        if test -x fxgraphics ; then \
	            echo "** Fortran77 application can be linked with graphics library" ; \
	            $(RM) -f fxgraphics ; \
	        else \
	            echo "** Fortran77 application CANNOT be linked with graphics library" ; \
	        fi ; \
	    fi ; \
	    echo ; \
	    $(MAKE) clean ; \
	  fi

clean:
	@-$(RM) -f work.pc work.pcl
	@-$(RM) -f *.o *~ PI* $(EXECS) \
	           cpilog cpi_anim cpi_trace cpi_log fpilog cxgraphics fxgraphics

distclean: clean
	@-$(RM) -f Makefile
