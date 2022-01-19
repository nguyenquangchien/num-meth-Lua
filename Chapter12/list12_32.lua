-- File list12_32.lua --
-- Example of capacitor with two dielectrics

require"pde2bv"

feqs = {
	function(x,y,uxx,ux,u,uy,uyy,i,j)
		if j==nt then ext = fyext*uy else ext = 0.0 end
		return uxx + uyy + ext
	end,
	function(x,u,uy,ux,i) return u end,
	function(x,u,uy,ux,i) return u-Vm end,
	function(y,u,ux,uy,j) return ux end,
	function(y,u,ux,uy,j) return ux end
}

ep1, ep2 = 3, 12	
Vm, xmax, ymax = 1.0, 1.0, 1.0
Nx,Ny = 80,80
nx,ny = Nx+1,Ny+1; nt = math.ceil(ny/2)
x,y = setxy({0,xmax},{0,ymax},Nx,Ny)
u = Spmat.new(nx,-ny) 
for j = 1,ny do	-- Set zero initial values
	u[j] = {}; for i = 1,nx do u[j][i] = 0 end
end
fyext = 2*(ep2-ep1)/((ep2+ep1)*(y[nt]-y[nt-1]))

SPM,COE,SOR,ADI = 1, 2, 3, 4 -- 4 solution types 
getfenv(pde2bvsor).nprint=1;getfenv(pde2bvcoe).nprint=1
u,errm = pde2bv(feqs,x,y,u,SOR) -- Replace SOR as desired
print('Dielectric change occurs at',nt)
print(' At end maximum correction =',errm)
ut = reversexy(u)
splot(ut); cplot(ut); splot('list12_32a.emf',ut)
write_data('list12_32a.dat',ut)
		
