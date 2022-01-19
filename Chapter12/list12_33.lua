-- File list12_33.lua --
-- Example of square corner resistor

require"pde2bv"

feqs = {
	function(x,y,uxx,ux,u,uy,uyy,i,j)
		ext = 0
		if i==nxmid then 
			if j<=nymid then ext = fxext*ux end
		end
		if j==nymid then
			if i<=nxmid then ext = ext + fyext*uy end
		end
		return uxx + uyy + ext
	end,
	function(x,u,uy,ux,i) 
		if i>nxmid-1 then return u
		else return uy end
	end,
	function(x,u,uy,ux,i) return uy end,
	function(y,u,ux,uy,j) 
		if j>nymid-1 then return u-Vm
		else return ux end
	end,
	function(y,u,ux,uy,j) return ux end
}
Vm = 1.0
Nx,Ny = 80,80
nx,ny = Nx+1,Ny+1
nxmid,nymid = math.ceil(nx/2),math.ceil(ny/2)

x,y = setxy({0,1},{0,1},Nx,Ny)
u = Spmat.new(nx,-ny) 
for j = 1,ny do	-- Set zero initial values
	u[j] = {}; for i = 1,nx do u[j][i] = 0 end
end
fxext = 4/(x[nxmid+1]-x[nxmid-1])
fyext = 4/(y[nymid+1]-y[nymid-1])

SPM,COE,SOR,ADI = 1, 2, 3, 4 -- 4 solution types 
getfenv(pde2bvsor).nprint=1;getfenv(pde2bvcoe).nprint=1
u,errm = pde2bv(feqs,x,y,u,COE) -- Replace COE as desired

print(' At end maximum correction =',errm)
ut = reversexy(u); splot(ut); cplot(ut)
ux,uy = grad(x,y,u)
sum = 0.5*(uy[1][nxmid] + uy[1][nx])
for i=nxmid+1,nx-1 do sum = sum + uy[1][i] end
sum = sum*(x[2]-x[1])
print('Eff squares =',1/sum)
splot('list12_33.emf',ut); cplot('list12_33a.emf',ut)
write_data('list12_33.dat',ut)
