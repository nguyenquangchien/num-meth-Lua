-- /* File list12_16.lua */
-- Program for transient solution of nonlinear transmission line

require"odefd"; require"pdeivbv"; require"intp"
getfenv(pdeivbv).nprint = 1
getfenv(ode2bvfd).umin = {5.e-4,5.e-4,5.e-4}
--getfenv(ode2bvfd).nprint = 1; getfenv(pdebivbv).nprint = 1

tt = os.time()
-- Model equations to be solved
us = 5.0 -- peak value of initial triangular wave

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

x = {}; nx = 1000; L = 20 -- Define x values 
for i=1,nx+1 do x[i] = L*(i-1)/nx end

u = {{},{}} -- Set initial wave values, 
for i=1,nx+1 do -- Set initial triangular wave at x=2
	if x[i]>3 then u[1][i] = 0
	elseif x[i]>2 then u[1][i] = us*(3-x[i])
	elseif x[i]>1 then u[1][i] = us*(x[i]-1)
	else u[1][i] = 0 end
end
for i=1,nx+1 do u[2][i] = intp(x,u[1],x[i],1) end

tvals = {0,2.5,25,20,{0,.25*L,.5*L,.75*L,L}} 
s,st = pdeivbvt({eq,efl,efr},tvals,x,u,up,upp)
write_data(2,'list12_16.dat',s) -- Save same number of points as before
write_data('list12_16t.dat',st)	
print('time taken for calculation = ',os.time()-tt)

