-- /* File list4_10.lua -- amplifier */
-- Programs to solve for the frequency response of a single stage amplifier

require"Complex"; require"nsolv" -- load complex,nsolv support
getfenv(nsolv).linear=1 -- Linear equations only need one iteration
-- Circuit Parameters
Rs,Rp,Re,Ru,Ro,Rc,Rl = 100, 20000, 200, 200000, 100000, 4000, 10000
Cs,Ce,Cc,Cp,Cu,Cce = 1.e-6, 10e-6, 2.e-6, 5.e-12, 1.e-12, 1.5e-12
gm = .01; vs = 1
-- Frequency factors
fact = 10^0.1-- 10 Points per decade in frequency
f = 1.0 -- Begin at frequency of 1Hz
jw = j*2*math.pi*f -- note jw is a variable name

eqs = function(y,v) -- Form of Node#1, Node#1, -- = Y[i] => Forced to be zero at solution values
	y[1] = (v[1]-vs)/(Rs+1/(jw*Cs)) + (v[1]-v[2])*(1/Rp+jw*Cp) + (v[1]-v[3])*(1/Ru+jw*Cu)
	y[2] = (v[2]-v[1])*(1/Rp+jw*Cp) + v[2]*(1/Re+jw*Ce) - gm*(v[1]-v[2]) + (v[2]-v[3])*(1/Ro+jw*Cce)
	y[3] = (v[3]-v[1])*(1/Ru+jw*Cu) + gm*(v[1]-v[2]) + (v[3]-v[2])*(1/Ro+jw*Cce) + (v[3]-v[4])*jw*Cc + v[3]/Rc
	y[4] = (v[4]-v[3])*jw*Cc + v[4]/Rl
end

sol = {{},{},{},{},{}} -- Table of 5 empty tables, filled later
local t1=os.clock() -- To check execution speed
nmax = 91; v = {0,0,0,0} -- initial guesses -- values not critical
for i=1,nmax do
	nsolv(eqs,v)
	v4 = v[4]; v4m = Complex.abs(v4)
	sol[1][i] = f
	sol[2][i] = math.log10(f)
	sol[3][i] = v4m
	sol[4][i] = Complex.angle(-v4)*180/math.pi - 180.
	sol[5][i] = 20*math.log10(v4m)
	f,jw = f*fact, jw*fact
end
print("time taken by nsolv = ",os.clock()-t1) -- Print time taken
write_data('list4_10.dat', unpack(sol)) -- Same as sol[1],sol[2],sol[3],sol[4]
plot(sol[2],sol[5]) -- Magnitude Bode plot
plot(sol[2],sol[4]) -- Angle Bode plot


