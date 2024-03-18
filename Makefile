FC=mpiifort
FOPTS=-O2 -w 
OBJS =  mocc.f input.f linita.f asymp.f ylog.f kmatrix.f diffcs.f odiffcs.f ediffcs.f outdat.f dpot.f potval.f math.f gauleg.f legendre.f
mocc: $(OBJS)
	$(FC) $(OBJS) -o $@ $(FOPTS)
	
