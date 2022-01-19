-- File list12_28.lua --
-- Example of BV problem in 2 dimensions  using pde2bv() function

require"pde2bv"
Vm = 5.0 -- Define equations to be solved
pi = math.pi;pi2 = 2*pi^2; sin=math.sin
xmax,ymax = 1,1
feqs = { -- Table of functions
	function(x,y,uxx,ux,u,uy,uyy,i,j) -- General point
		return uxx + uyy +Vm*pi2*(sin(pi*x/xmax)*sin(pi*y/ymax))^20
	end,
	function(x,u,uy,ux,i) return u end, -- Bottom boundary
	function(x,u,uy,ux,i) return u end, -- Top boundary
	function(y,u,ux,uy,j) return u end, -- Right boundary
	function(y,u,ux,uy,j) return u end  -- Left boundary
}  -- End general point and boundary values

Nx,Ny = 100,100
nx,ny = Nx+1,Ny+1
x,y = setxy({0,xmax},{0,ymax},Nx,Ny)
u = Spmat.new(nx,-ny) 
for j = 1,ny do	-- Set zero initial values
	u[j] = {}; for i = 1,nx do u[j][i] = 0 end
end

SPM,COE,SOR,ADI = 1, 2, 3, 4 -- 4 solution types 
getfenv(pde2bv).nprint=1; getfenv(pde2bvspm).nprint=1
t1 = os.clock()
u,errm = pde2bv(feqs,x,y,u,SPM) -- Replace ADI as desired
print('time =',os.clock()-t1)
print(' At end maximum correction =',errm)
splot(u); cplot(u)
splot('list12_28.emf',u)
write_data('list12_28.dat',u)
