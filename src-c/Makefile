LIBNAME=acepack
LD=ld

OBJS=	rlsmo.o \
	ace.o \
	avas.o 

acepack.so:	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS) 
	echo "If acepack does not load properly add"
	echo " -lf2c to the Makefile line beginning"
	echo "  @$(LD) "

clean:
	@rm -f *.o *.so

