-- /* File list12_3.lua */
-- Programs to integrate diffusion equation in one spatial variable

require"pdeivbv"
getfenv(odefd).nprint = 1; getfenv(pdeivbv).nprint = 1
-- Model equations to be solved
L,D,U0 = 1.0,1.0,1.0; u = {}

eq = function(fu,x,t,u,up,upp,ut) -- Equations with time and spatial derivatives
	fu[1] = ut[1] - D*upp[1] -- Diffusion equation
end
efl = function(fu,u,up) -- Left boundary fixed at 1.0
	fu[1] = u[1] - U0
end
efr = function(fu,u,up) -- Right boundary, fixed at 0.0
	fu[1] = u[1] - 0.0  
end

x = {}; Nx = 400; tmax = 0.01
for i=1,Nx+1 do x[i],u[i] = (i-1)*L/Nx, 0.0 end -- Set up x and u arrays

sol = pdeivbv({eq,efl,efr},{0,tmax,5,5},x,u)
plot(sol[1],sol[2],sol[3],sol[#sol])
write_data('list12_3.dat',sol)
