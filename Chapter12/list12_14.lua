-- /* File list12_14.lua */
-- Program for transient solution of transmission line

require"odefd"; require"pdeivbv"
getfenv(pdeivbv).nprint = 1
getfenv(ode2bvfd).umin = {5.e-4,5.e-4,5.e-4}
--getfenv(ode2bvfd).nprint = 1; getfenv(pdebivbv).nprint = 1

tt = os.time()
-- Model equations to be solved
Zo = 50
Vt, Is = .026, 1.e-12 -- Thermal voltage, Sat current
Rr = 0 -- Series resistance factor
Rg = 0 -- Parallel conductance factor
Rs = 100 -- Source resistance
Rl = 100 -- Load resistance
vs = 5.0 -- Step value of source voltage
w = 2*math.pi; Rrg, Rrpg = Rr*Rg, Rr + Rg

eq = function(fv,x,t,v,vp,vpp,vt,vtt) -- Equations with time and spatial derivatives
	fv[1] = vpp[1] - Rrg*v[1] - Rrpg*vt[1] - vtt[1] -- Second order equations
	fv[2] = vpp[2] - Rrg*v[2] - Rrpg*vt[2] - vtt[2]
	end

efl = function(fv,v,vp,t,vt) -- Using mixed boundary condition
	fv[1] = vs*math.sin(w*t) - v[1] - Rs*v[2]
	fv[2] = vp[1] + Zo*(Rr*v[2] + vt[2])
end
efr = function(fv,v,vp,t,vt) -- Using mixed boundary condition
	fv[1] = v[2] - Is*(math.exp(v[1]/Vt) - 1) -- Diode at load
	fv[2] = Zo*vp[2] + Rg*v[1] + vt[1]
end

x = {}; nx = 500; L = 1 -- Define x values 
for i=1,nx+1 do x[i] = L*(i-1)/nx end

v = {{},{}} -- Set initial voltage values, V = 0, I = 0 for all x values
for i=1,nx+1 do
	v[1][i] = 0; v[2][i] = 0 -- v[2] is current variable
end
tvals = {0,10,100,40,{0,.25,.5,.75,1}} 
s,st = pdeivbvt({eq,efl,efr},tvals,x,v,vp,vpp)
write_data(2,'list12_14.dat',s) -- Save same number of points as before
write_data('list12_14t.dat',st)	
print('time taken for calculation = ',os.time()-tt); tt = os.time()

