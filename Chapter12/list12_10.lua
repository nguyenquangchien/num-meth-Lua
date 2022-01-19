-- /* list12_10.lua */
-- Solution of nonlinear Chemical Engineering IV-BV problem
require"pdeivbv"
getfenv(pdeivbv).nprint = 1

M,EdR,b,Q,ce,Te = 2, 22000, .5e8, 1000, 0.07, 1250
Ti = 1270 -- Initial temperature
L = 48; Nx = 200; Tm = 60; Nt = 6 -- spatial, time parameters

f = function(eqs,x,t,u,up,upp,ut) -- Differntial equation
	eqs[1] = upp[1]/M - up[1] - b*u[1]^2*math.exp(-EdR/u[2]) - ut[1]
	eqs[2] = upp[2]/M - up[2] +Q*b*u[1]^2*math.exp(-EdR/u[2]) - ut[2]
end
fl = function(eqs,u,up) -- Left boundary conditions
	eqs[1] = u[1] - up[1]/M - ce
	eqs[2] = u[2] -up[2]/M - Te
end
fr = function(eqs,u,up)
	eqs[1] = up[1] -- Zero slopes
	eqs[2] = up[2]
end

tvals = {0,Tm,Nt}
x,u = {},{{},{}} -- Set initial values
for i=1,Nx+1 do x[i],u[1][i],u[2][i] = L*(i-1)/Nx,0,1270 end
s = pdeivbv({f,fl,fr},tvals,x,u) -- Solve equations
plot(s[1],s[#s-1]); plot(s[1],s[#s])
write_data('list12_10.dat',s)

