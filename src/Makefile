LIBNAME=acepack
LD=ld

OBJS=	rlsmo.o \
	ace.o \
	avas.o 

acepack.so:	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so

