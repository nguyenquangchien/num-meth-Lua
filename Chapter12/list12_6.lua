-- /* File list12_6.lua */
-- Programs to integrate diffusion equation in one spatial variable

require"pdeivbv"; require'odebvfd'
require"elemfunc"; erfc = elemfunc.erfc
getfenv(pdeivbv).nprint = 1; getfenv(pdeivbvqs).nprint=1
getfenv(ode2bvfd).NMAX = 1 -- Only one iteration required
-- Model equations to be solved
L,D = 1,1; u = {}

eq = function(fu,x,t,u,up,upp,ut) -- Equations with time and spatial derivatives
	fu[1] = ut[1] - D*upp[1] -- Diffusion equation
end
efl = function(fu,u,up) -- Left boundary fixed at 1.0
	fu[1] = u[1] - 1.0
end
efr = function(fu,u,up) -- Right boundary, fixed at 0.0
	fu[1] = u[1] - 0.0  
end

Nx = 400
--Nx = 200 -- Half spatial steps
x = xlg(0,1,1.e-6,Nx) -- Set up x and u arrays
for i=1,Nx+1 do u[i] = 0.0 end 

sol = pdeivbvqs({eq,efl,efr},{0,{1.e-6,1.e-2},1,20},x,u)
--sol = pdeivbvqs({eq,efl,efr},{0,{1.e-6,1.e-2},1,10},x,u) -- Half time steps
u1,u2 = {},{} -- Theoretical values
t = 1.e-6
for i=1,Nx+1 do u1[i] = erfc(x[i]/(2*math.sqrt(t))) end 
t = 1.e-2
for i=1,Nx+1 do u2[i] = erfc(x[i]/(2*math.sqrt(t))) end 
emax = 0.0
for i=1,Nx+1 do
	 ex = math.abs(sol[3][i] - u1[i])
	 if ex>emax then emax = ex end
end
print('Max error at t = 1.e-6 is ',emax)
emax = 0.0
for i=1,Nx+1 do
	 ex = math.abs(sol[7][i] - u2[i])
	 if ex>emax then emax = ex end
end
print('Max error at t = 1.e-2 is ',emax)
plot(sol[1],sol[2],sol[3],sol[#sol])
write_data('list12_6.dat',sol,u1,u2)
