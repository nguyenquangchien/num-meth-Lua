-- /* File list12_7.lua */
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

Nx = 400; Nt = 20
x1 = xlg(0,1,1.e-6,Nx) -- Set up x and u arrays
for i=1,Nx+1 do u[i] = 0.0 end 
sol1 = pdeivbvqs({eq,efl,efr},{0,{1.e-6,1.e-2},1,Nt},x1,u)

Nx = 200; Nt = 10 -- Half spatial steps, half time steps
x2 = xlg(0,1,1.e-6,Nx) -- Set up x and u arrays
for i=1,Nx+1 do u[i] = 0.0 end 
sol2 = pdeivbvqs({eq,efl,efr},{0,{1.e-6,1.e-2},1,Nt},x2,u)

err1 = odeerror({x1,sol1[3]},{x2,sol2[3]}) -- t = 1.e-6
err2 = odeerror({x1,sol1[7]},{x2,sol2[7]}) -- t = 1.e-2
write_data('list12_7.dat',err1,err2)
