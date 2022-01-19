-- File list13_18.lua -- Capacitor with nonlinear dielectric

require'pde2fe' 

nts = {readnts('list13_11')}

Vm = 1.0
ep1,gr = 3,3
feqs = {
	function(x,y,uxx,ux,u,uy,uyy,ntr,mtn)
		ep = ep1*(1 + gr*y)^2*(1+gr*u^2)
		return ep*(uxx + uyy) -- OK 
	end,
	function(nd,u,un,nbs,kb)
		if nbs==1 then return u -- Bottom boundary
		elseif nbs==3 then return u - Vm -- Top boundary
		else  return un end -- Sides 
	end 
}	

getfenv(pde2fe).nprint = 1 -- Observe progress

u = pde2fe(feqs,nts)

x,y = {},{}; NT = 21
for j=1,NT do y[j] = (j-1)/(NT-1); x[j] = 2*y[j] end
sol = toxysol(u,x,y,nts[1],nts[2])

write_data('list13_18.dat',sol)
splot(sol)
