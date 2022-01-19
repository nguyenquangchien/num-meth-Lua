-- File list13_13.lua -- FE analysis of Square Corner Resistor

require'pde2fe' 

nts = {readnts('list13_13')}
Un = 1 -- Normal derivative on left boundary

feqs = {
	function(x,y,uxx,ux,u,uy,uyy)
		return uxx + uyy
	end,
	function(nd,u,un,nbs)
		if nbs==1 then return u -- Bottom boundary
		elseif nbs==7 then return u -- Bottom end points
		elseif nbs==4 then return un - Un -- Left boundary
		else  return un end -- Sides 
	end 
}	

getfenv(spgauss).nprint = 1 -- Just to see progress
getfenv(pde2fe).linear = 1 -- Saves an extra Newton iteration
u = pde2fe(feqs,nts)

pt = {0,1.5} -- Point at center of left boundary
print('Effective squares of corner = ',intptri(u,nts[1],nts[2],pt)-2)
x,y = {},{}; NT1 = 21; NT2 = 41 -- Arrays for plotting
for j=1,NT1 do x[j] = 1+(j-1)/(NT1-1) end 
for j=1,NT2 do y[j] = 2*(j-1)/(NT2-1) end
sol1 = toxysol(u,x,y,nts[1],nts[2]); splot(sol1)
write_data('list13_13a.dat',sol1) -- Save for 1<x<2; 0<y<2
for j=1,NT1 do x[j] = (j-1)/(NT1-1); y[j] = 1+x[j] end
sol2 = toxysol(u,x,y,nts[1],nts[2]); splot(sol2)
write_data('list13_13b.dat',sol2) -- Save for 0<x<1; 1<y<2