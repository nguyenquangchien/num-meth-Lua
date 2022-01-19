-- /* File list12_17.lua */
-- Program for transient solution of nonlinear transmission line

require"odefd"; require"pdeivbv"
require"intp"
getfenv(pdeivbv).nprint = 1
getfenv(ode2bvfd).umin = {5.e-4,5.e-4,5.e-4}
--getfenv(ode2bvfd).nprint = 1; getfenv(pdebivbv).nprint = 1

tt = os.time()
-- Model equations to be solved
cosh = function(x) return (math.exp(x) + math.exp(-x))*.5 end
ui = function(x,a,xo) return a/(cosh(0.5*math.sqrt(2*a)*(x-xo))^2) end

eq = function(fu,x,t,u,up,upp,ut) -- Equations with time and spatial derivatives
	fu[1] = up[1] - u[2]
	fu[2] = upp[2] + 6*u[1]*u[2] + ut[1]
end
efl = function(fu,u,up,t,ut) -- Left boundary condition
	fu[1], fu[2] = u[1], u[2]
end
efr = function(fu,u,up,t,ut) -- Right boundary condition
	fu[1], fu[2] = u[1], u[2]
end

x = {}; nx = 1000; L = 40 -- Define x values 
for i=1,nx+1 do x[i] = L*(i-1)/nx end

u = {{},{}} -- Set initial two solitons at x=4, x=20 
for i=1,nx+1 do u[1][i] = ui(x[i],4,4) end
for i=1,nx+1 do u[1][i] = u[1][i] + ui(x[i],1,20) end
for i=1,nx+1 do u[2][i] = intp(x,u[1],x[i],1) end

tvals = {0,5,50,10,{0,.25*L,.5*L,.75*L,L}} 
s,st = pdeivbvt({eq,efl,efr},tvals,x,u)
write_data(2,'list12_17.dat',s) -- Save same number of points as before
write_data('list12_17t.dat',st)	
print('time taken for calculation = ',os.time()-tt); tt = os.time()

