-- File list13_17.lua -- Deflection of circular plate fixed at edges.

require'pde2fe' 
L = 1.0 -- Radius of circular plate
nts ={ readnts('list13_17')}
nds,tri = nts[1], nts[2] -- Node data

E,dz,sig,q = 2e11, 0.01, 0.3, 3.36e4 -- Parameters
D = E*dz^3/(12*(1-sig^2)); print('D = ',D)

feqs = {
	function(x,y,uxx,ux,u,uy,uyy,ntr,mtn)
		return uxx + uyy  - q/D
	end,
	function(nd,u,un,nbs)
		return u
	end
}

getfenv(spgauss).nprint = 1; getfenv(pde2fe).nprint = 1
getfenv(pde2fe).linear = 1 -- Observe progress
uv = pde2fe(feqs,nts) -- Solve equations, using SPM

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
z = pde2fe(feqz,nts) -- Try other solution methods?

x,y = {},{}; NT = 41 -- x,y grid for plotting results
for j=1,NT do x[j] = 2*L*(j-1)/(NT-1)-L; y[j] = x[j] end
sol = toxysol(uv,x,y,nds,tri); solx = toxysol(z,x,y,nds,tri)
write_data('list13_17.dat',sol,solx)
print('Maximum deflection =',intptri(z,nds,tri,{0,0}))
splot(sol); splot(solx)