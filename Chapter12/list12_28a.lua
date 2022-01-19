-- File list12_28.lua --
-- Example of BV problem in 2 dimensions  using pde2bv() function

require"pde2bv"
Vm = 5.0 -- Define equations to be solved
pi = math.pi;pi2 = 2*pi^2; sin=math.sin
feqs = { -- Table of functions
	function(x,y,uxx,ux,u,uy,uyy,i,j) -- General point
		return uxx + uyy + Vm*pi2*(sin(pi*x)*sin(pi*y))^20
	end,
	function(x,u,uy,ux,i) return u end, -- Bottom boundary
	function(x,u,uy,ux,i) return u end, -- Top boundary
	function(y,u,ux,uy,j) return u end, -- Right boundary
	function(y,u,ux,uy,j) return u end  -- Left boundary
}  -- End general point and boundary values

Nx,Ny = 100,100
nx,ny = Nx+1,Ny+1
x,y = setxy({0,1},{0,1},Nx,Ny)
u = Spmat.new(nx,-ny) 
for j = 1,ny do	-- Set zero initial values
	u[j] = {}; for i = 1,nx do u[j][i] = 0 end
end

SPM,COE,SOR,ADI = 1, 2, 3, 4 -- 4 solution types 
getfenv(pde2bvadi).nprint=1; getfenv(pde2bvcoe).nprint=1
u,errm = pde2bv(feqs,x,y,u,COE) -- Replace ADI as desired

print(' At end maximum correction =',errm)
splot(u); cplot(u)
splot('list12_28.emf',u)
write_data('list12_28.dat',u)
