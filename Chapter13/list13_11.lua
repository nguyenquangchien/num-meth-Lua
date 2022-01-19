-- File list13_11.lua -- 2 dielectric capacitor

require'pde2fe' 

nts = {readnts('list13_11')}

Vm = 1.0
ep1,ep2 = 3, 12
feqs = {
	function(x,y,uxx,ux,u,uy,uyy,ntr,mtn)
		if mtn==1 then ep = ep1 
		else ep = ep2 end
		return ep*(uxx + uyy)
	end,
	function(nd,u,un,nbs,kb)
		if nbs==1 then return u -- Bottom boundary
		elseif nbs==3 then return u - Vm -- Top boundary
		else  return un end -- Sides 
	end 
}	

getfenv(spgauss).nprint = 1
getfenv(pde2fe).nprint = 1
getfenv(pde2fe).linear = 1

u = pde2fe(feqs,nts)

x,y = {},{}; NT = 21
for j=1,NT do y[j] = (j-1)/(NT-1); x[j] = 2*y[j] end
sol = toxysol(u,x,y,nts[1],nts[2])

write_data('list13_11.dat',sol)
splot(sol)
