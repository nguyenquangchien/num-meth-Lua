-- /* File list12_8.lua */
-- Programs to integrate nonlinear diffusion equation in one spatial variable

require"pdeivbv"; require'odebvfd'
getfenv(pdeivbv).nprint = 1; getfenv(pdeivbvqs).nprint=1
--getfenv(ode2bvfd).nprint=1 -- See detailed convergence

-- Model equations to be solved
D00,D10,E0,E1 = 0.05, 0.95, 3.5, 3.5 -- Diffusion coeff parameters
T = 1000+273 -- temperature
D0sav = D00*math.exp(-E0/(0.026*T/300))
D1sav = D10*math.exp(-E1/(0.026*T/300))
ni = 7.14e18; L = 1.e-4 ; ul = 5e20; ur = 0; u = {}
print('D0, D1 =',D0sav,D1sav)

eq = function(fu,x,t,u,up,upp,ut) -- Equations with time and spatial derivatives
	D = D0 + D1*u[1]/ni
	fu[1] = ut[1] - D*upp[1] - D1*up[1]^2/ni
end
efl = function(fu,u,up) -- Left boundary 
	fu[1] = u[1] - ul 
end
efr = function(fu,u,up) -- Right boundary
	fu[1] = u[1] - ur  
end

Nx = 200; Nt = 20 -- spatial steps, time steps
x = xlg(0,L,L*1.e-4,Nx) -- Set up x and u arrays
for i=1,Nx+1 do u[i] = 0.0 end; u[1] = ul 

D1=0; D0 = D0sav+D1sav*(ul/ni) -- Get initial trial solution
pdeivbv({eq,efl,efr},{0,.001,1},x,u) -- Constant diffusion coeff
D1=D1sav; D0 = D0sav -- Now real diffusion coeff
sol = pdeivbvqs({eq,efl,efr},{0,{.01,1.e3},1,Nt},x,u)
plot(sol[1],sol[#sol])
write_data('list12_8.dat',sol)
