-- /* list3_7.lua */ -- Newton's method solution for diode circuit

require "newton"

Vs,R,Is,vt = 15.00,1.e3,1.e-14,.026 -- Circuit and device parameters
xa,vo,vsa = {},{},{} -- Define empty arrays to hold calculated values

-- Rectifier circuit equation to be solved
function fv(v) return vs - v - R*Is*(math.exp(v/vt) - 1) end 

nmx = 0
for i=1,201 do -- Select 200 intervals for calculation
	x = (i-1)*4*math.pi/(200) -- Angle value
	vs = Vs*math.sin(x) -- source voltage
	vd,n,err = newton(fv, math.min(vs,0.8)) -- diode voltage by Newton's method
	nmx = math.max(nmx,n) -- Check on convergence
	xa[i], vo[i], vsa[i] = x, vs-vd, vs -- Save arrays
end
print('Maximum number of iterations =',nmx)
write_data('list3_7a.dat',xa,vsa,vo)
plot(xa,vsa,vo)