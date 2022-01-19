-- /* File list12_5.lua */
-- Programs to integrate diffusion equation in one spatial variable

require"pdeivbv"; require'odebvfd'
getfenv(odefd).nprint = 1; 
getfenv(pdeivbv).nprint = 1
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

Nx = 400; x = xlg(0,1,1.e-6,Nx) -- Set up x and u arrays
for i=1,Nx+1 do u[i] = 0.0 end 

sol = pdeivbvqs({eq,efl,efr},{0,{1.e-6,1.e-1},1,20},x,u)
plot(sol[1],sol[2],sol[3],sol[#sol])
write_data('list12_5.dat',sol)
