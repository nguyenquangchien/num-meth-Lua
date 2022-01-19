-- File list13_6.lua --

require'pde2fe' 

nts = {readnts('list13_2')}

pi = math.pi; pi2 = 2*pi^2; sin = math.sin
feqs = {
	function(x,y,uxx,ux,u,uy,uyy,ntr,mtn)
		return uxx + uyy  + pi2*(sin(pi*x)*sin(pi*y))
	end,
	function(nd,u,un,nbs)
		return u
	end
}	

getfenv(spgauss).nprint = 1; getfenv(pde2fe).nprint=1
getfenv(pde2fe).linear = 1

u = pde2fe(feqs,nts)
--u = pde2fe(feqs,nts,_,3) -- Try approximate methods

x,y = {},{}; NT = 41
for j=1,NT do x[j] = (j-1)/(NT-1); y[j] = x[j] end
sol,solxy = toxysol(u,x,y,nts)

write_data('list13_6a.dat',sol); write_data('list13_6b.dat',solxy)
splot(solxy); cplot(solxy)