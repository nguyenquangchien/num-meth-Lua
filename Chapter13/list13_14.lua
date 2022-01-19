-- File list13_14.lua -- Thermal heating of square corner resistor

require'pde2fe' 

Um,K = 10.0,6.9; L = 1.e-2 -- L in cm
nts = {readnts('list13_13',L)} -- Read spatial grid; Scale to L size
getfenv(spgauss).nprint = 1 -- just to see progress
getfenv(pde2fe).linear = 1 -- saves an extra Newton iteration

fels = { -- Equations for electric field
	function(x,y,uxx,ux,u,uy,uyy,ntr,mtn)
		return uxx + uyy 
	end,
	function(nd,u,un,nbs,kb)
		if nbs==1 then return u -- Bottom boundary
		elseif nbs==7 then return u -- Bottom end points
		elseif nbs==4 then return u - Um -- Left boundary
		elseif nbs==9 then return u - Um
		else  return un end -- All other sides 
	end 
}
uel = pde2fe(fels,nts) -- Solve for electric potential

efdx,efdy = derivtri(nts[2],nts[1],uel) -- Get electric field
esq = {} -- Square of electric field
for i=1,#efdx do esq[i] = efdx[i]^2+efdy[i]^2 end

feqs = { -- table of PDE and boundary functions for temperature
	function(x,y,uxx,ux,u,uy,uyy,ntr,mtn)
		return uxx + uyy + K*esq[ntr]
	end,
	function(nd,u,un,nbs,kb)
		if nbs==1 then return u -- Bottom boundary
		elseif nbs==7 then return u -- Bottom end points
		elseif nbs==4 then return u -- Left boundary
		elseif nbs==9 then return u -- Left end points
		else  return un end -- All other sides 
	end 
}	
T = pde2fe(feqs,nts) -- Solve equations for temperature

x,y = {},{}; NT1 = 21; NT2 = 41 -- Rectangular array for plotting
for j=1,NT1 do x[j] = L*(1+(j-1)/(NT1-1)) end
for j=1,NT2 do y[j] = 2*L*(j-1)/(NT2-1) end
sol1 = toxysol(T,x,y,nts[1],nts[2]); splot(sol1)
for j=1,NT1 do x[j] = L*(j-1)/(NT1-1); y[j] = L+x[j] end
sol2 = toxysol(T,x,y,nts[1],nts[2]); splot(sol2)
Tm = 0; for i=1,#T do if T[i]>Tm then Tm,im = T[i],i end end
print('Maximum T =',Tm,' at x,y =',nts[1][im][1],nts[1][im][2])
write_data('list13_14a.dat',sol1)
write_data('list13_14b.dat',sol2)

