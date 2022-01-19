-- File list12_31.lua --
-- Example of BV problem for capacitor with offset bottom plate

require"pde2bv"
local xmax,ymax
Vm = 1
feqs = { -- Function and boundary conditions
	function(x,y,uxx,ux,u,uy,uyy,i,j)
		if x>x1 then
			if y<=y1 then return u else return uxx+uyy end
		else	return uxx + uyy end
	end,
	function(x,u,uy,ux,i) return u end,
	function(x,u,uy,ux,i) return u - Vm end,
	function(y,u,ux,uy,j) return ux end,
	function(y,u,ux,uy,j) 
		if y<=y1 then return u else return ux end
	end
}

x,y = {},{}; Nx,Ny = 40,40
nx,ny = Nx+1,Ny+1
xmax,ymax = 5, 1
x1,y1 = xmax/2, ymax/2
for i=1,nx do x[i] = xmax*(i-1)/Nx end
for j=1,ny do y[j] = ymax*(j-1)/Ny end
u = Spmat.new(nx,-ny)
for j = 1,ny do	
	u[j] = {}; for i = 1,nx do u[j][i] = 0 end
end

SPM,COE,SOR,ADI = 1, 2, 3, 4 -- 4 solution types 
getfenv(pde2bvsor).nprint=1
u,errm = pde2bv(feqs,x,y,u,SOR) -- Replace SOR as desired
ur = reversexy(u); splot(ur); cplot(ur)
ux,uy = grad(x,y,u)
sum = 0.5*(uy[ny][1] + uy[ny][nx])
for i=2,nx-1 do sum = sum + uy[ny][i] end
sum = sum*(x[2]-x[1])
print('C/eps*W =',sum)
write_data("list12_31.dat",ur); cplot("list12_31.emf",ur)


