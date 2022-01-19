-- /* File list4_11.lua  */
-- Programs #2 for the frequency response of an amplifier

require"Complex"; require"gauss" -- load complex variable support
--getfenv(nsolv).full = true -- Try gauss() instead of spgauss
-- Circuit Parameters
Rs,Rp,Re,Ru,Ro,Rc,Rl = 100, 20000, 200, 200000, 100000, 4000, 10000
Cs,Ce,Cc,Cp,Cu,Cce = 1.e-6, 10e-6, 2.e-6, 5.e-12, 1.e-12, 1.5e-12
gm = .01; vs = 1
-- Frequency factors
fact = 10^0.1-- 10 Points per decade in frequency
f = 1.0 -- Begin at frequency of 1Hz
jw = j*2*math.pi*f

sol = {{},{},{},{},{}} -- Table of 5 empty tables
local c=os.clock()
nmax = 91; v = {0,0,0,0} -- initial guesses -- values not critical
local t1=os.clock()
-- Second approach -- Set up matrix equations 
for i=1,nmax do -- gmat is conductance matrix of coefficients
	gmat = {{(1/(Rs+1/(jw*Cs))+1/Rp+jw*Cp+1/Ru+jw*Cu), (-1/Rp-jw*Cp), (-1/Ru-jw*Cu),0},
		{(-1/Rp-jw*Cp-gm),(1/Rp+jw*Cp+1/Re+jw*Ce+gm+1/Ro+jw*Cce),(-1/Ro-jw*Cce),0},
		{(-1/Ru-jw*Cu+gm),(-gm-1/Ro-jw*Cce),(1/Ru+jw*Cu+1/Ro+jw*Cce+jw*Cc+1/Rc),(-jw*Cc)},
		{0,0,(-jw*Cc),(1/Rl+jw*Cc)}}
	b = {vs/(Rs+1/(jw*Cs)),0,0,0} -- Source terms -- Dimension of current
	gauss(gmat,b) -- Solutions for v's are returned in the b array 
	v4 = b[4]; v4m = Complex.abs(v4)
	sol[1][i] = f
	sol[2][i] = math.log10(f)
	sol[3][i] = v4m
	sol[4][i] = Complex.angle(-v4)*180/math.pi -180.
	sol[5][i] = 20*math.log10(v4m)
	f,jw = f*fact, jw*fact
end
print("time taken by gauss = ",os.clock()-t1)
write_data('list4_11.dat', unpack(sol))
plot(sol[2],sol[5]) -- Magnitude Bode plot
plot(sol[2],sol[4]) -- Angle Bode plot


