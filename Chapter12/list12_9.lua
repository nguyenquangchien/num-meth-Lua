-- /* File list12_9.lua */
-- Programs to integrate nonlinear diffusion equation in one spatial variable

require'odebvfd'
-- Model equations to be solved
D00,D10,E0,E1 = 0.05, 0.95, 3.5, 3.5 -- Diffusion coeff parameters
T = 1000+273 -- temperature
D0 = D00*math.exp(-E0/(0.026*T/300))
D1 = D10*math.exp(-E1/(0.026*T/300))
ni = 7.14e18; ur = 0.0; x,u = {},{}
L = 43; ul = 1e21; ex = 'a' -- Change as desired
--L = 30; ul = 5e20; ex = 'b'
--L = 14; ul = 1e20; ex = 'c'
--L = 10; ul = 5e19; ex = 'd'
--L = 6; ul = 1e19; ex = 'e'
--L = 4; ul = 5e18; ex = 'f'
--L = 4; ul = 1e18; ex = 'g'
--L = 4; ul = 1e17; ex = 'h'
Dr = D1*ul/(D0*ni)
print(1/(2*math.sqrt(D0*1e3)),D0,D1)

eq = function(x,u,up,upp) -- Equations with time and spatial derivatives
	D = 1 + Dr*u
	return D*upp + (Dr*up + 2*x)*up
end

Nx = 10000 -- spatial steps, time steps
-- Set up x and u arrays
for i=1,Nx+1 do x[i] = (i-1)*L/Nx end 
u,nn,err1,err2 = ode1fd(eq,x,{1.0,0})

plot(u); write_data('list12_9'..ex..'.dat',u)
