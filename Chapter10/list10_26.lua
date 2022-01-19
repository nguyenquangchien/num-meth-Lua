--- /* File list10_26.lua */
-- Programs to solve AC voltage multiplier circuit
require"odeiv"

C1,C2,C3,C4,Rl = 100e-6, 100e-6, 100e-6, 100e-6, 100E3
Vm,w,tmax = 10.0, 2*math.pi*60, 0.7
Is,vth = 1.e-12, 0.026

id = function(v)
	return Is*(math.exp(v/vth) - 1)
end
f = function(eqs,t,v,vp)
	eqs[1] = C1*(vp[1]-vp[5]) + C2*(vp[1]-vp[3]) - id(-v[1]) + id(v[1]-v[2]) 
	eqs[2] = C3*vp[2] + C4*(vp[2]-vp[4]) + id(v[2]-v[3]) - id(v[1]-v[2])
	eqs[3] = C2*(vp[3]-vp[1]) - id(v[2]-v[3]) + id(v[3]-v[4])
	eqs[4] = C4*(vp[4]-vp[2]) - id(v[3]-v[4]) + v[4]/Rl
	eqs[5] = v[5] - Vm*math.sin(w*t)
end
s1,err = odeivse(f,{0,tmax},{0,0,0,0,0})

plot(s1[1],s1[5]); plot(err[1],err[5])
write_data("list10_26.dat",s1,err)

