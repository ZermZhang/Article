#求解近似精确解的Make文件
MGD9V	= ../comm/mgd9v.f
COEF	= ../comm/coef.f

lin2048 : linear2048.f $(MGD9V) $(COEF)
	f77 -o lin2048 linear2048.f $(MGD9V) $(COEF) -mcmodel=medium
	
lin4096 : linear4096.f $(MGD9V) $(COEF)
	f77 -o lin4096 linear4096.f $(MGD9V) $(COEF) -mcmodel=medium

run2048 : lin2048
	./lin2048
	
run4096 : lin4096
	./lin4096
	
	
clean :
	rm lin2048 *~
	rm lin4096 *~
