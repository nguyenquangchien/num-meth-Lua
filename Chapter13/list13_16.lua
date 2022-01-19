-- File list13_16.lua -- Deflection of square plate fixed at edges.

require'pde2fe' 
L = 2.0 -- Scale length to 2 by 2
nts ={ readnts('list13_15',L)} -- Square area with 2996 nodes
nds,tri = nts[1], nts[2] -- Node data

E,dz,sig,q = 2e11, 0.01, 0.3, 3.36e4
D = E*dz^3/(12*(1-sig^2)); print('D = ',D)

feqs = {
	function(x,y,uxx,ux,u,uy,uyy)
		return uxx + uyy  - q/D
	end,
	function(nd,u,un,nbs)
		return u
	end
}

getfenv(spgauss).nprint = 1; getfenv(pde2fe).nprint = 1
getfenv(pde2fe).linear = 1 -- Saves Newton iterations

uv = pde2fe(feqs,nts) -- Solve equations for uz

feqz = {
	function(x,y,uxx,ux,u,uy,uyy,ntr,mtn,nds,ltr)
		local n1,n2,n3 = nds[1],nds[2],nds[3]
		local l1,l2,l3 = ltr[1],ltr[2],ltr[3]
		return uxx + uyy  - (uv[n1]*l1+uv[n2]*l2+uv[n3]*l3)
	end,
	function(nd,u,un,nbs)
		return u
	end
}
z = pde2fe(feqz,nts) -- Solve for deflection

x,y = {},{}; NT = 41 -- x - y grid for plotting
for j=1,NT do x[j] = L*(j-1)/(NT-1); y[j] = x[j] end
sol = toxysol(uv,x,y,nds,tri); solx = toxysol(z,x,y,nds,tri)
write_data('list13_16.dat',sol,solx)
splot(sol); splot(solx)