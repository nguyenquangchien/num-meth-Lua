-- File list13_24.lua -- Vibration of circular membrane.

require'pde2fe' ; require'elemfunc'
L = 1.0 -- Radius of circular plate
nts ={ readnts('list13_17')}
nds,tri = nts[1], nts[2] -- Node data

C = 1; Csq = C^2-- Velocity -- change as desired
gm = 0.0 -- damping coefficient
--gm = .2 -- change as desired
x,y = {},{}; NT = 41
for j=1,NT do x[j] = 2*L*(j-1)/(NT-1)-L; y[j] = x[j] end

feqs = {
	function(x,y,t,uxx,ux,u,uy,uyy,ut,utt)
		return uxx + uyy  - utt/Csq - gm*ut/C -- Wave equation in 2D
	end,
	function(nd,t,u,un,ut,utt) -- Boundary value
		return u
	end
}

u = {}
for k=1,#nds do
	r = math.sqrt(nds[k][1]^2+nds[k][2]^2)
	u[k] = elemfunc.J0(5.52008*r) -- Use Bessel function 
end
SPM,COE,SOR = 1, 2, 3 -- 3 solution methods
getfenv(pde2fe1t).nprint=1
getfenv(pde2fe).linear = 1
tvals = {0,8,200,1,{{0,0},{.5,0},{-.5,0},{.5,.5},{-.5,-.5},{.75,0},{0,.75}}}
u,uxyt = pde2fe1t(feqs,tvals,nts,u,COE)

nsol = #u
for i=1,nsol,10 do -- Save 2D data files
	sfl = 'list13_24.'..(1+(i-1)/10)..'.dat'
	sol,solxy = toxysol(u[i],x,y,nts[1],nts[2])
	write_data(sfl,sol); splot(sol)
end
write_data('list13_24a.dat',uxyt)
