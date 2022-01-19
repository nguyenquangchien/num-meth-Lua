--- /* File list10_28.lua */
-- Programs to solve for transient response of amplifier
require"odeiv"

--Circuit Parameters
Rs,Rp,Re,Ru,Ro,Rc,Rl = 100, 20e3, 200, 200e3, 100e3, 4e3, 10e3
Cs,Ce,Cc,Cp,Cu,Cce = 1e-6, 10e-6, 2e-6, 5e-12, 1e-12, 1.5e-12
gm = 0.01
Vm,t1,t2 = -.1, 2.e-7, 4.e-7 -- Short pulse (200, 400 nsec)

vs = function(t)
	if t<t1 then return Vm else return 0 end
end
f = function(eqs,t,v,vp)
	eqs[1] = Cs*(vp[1]-vp[5])+Cu*(vp[1]-vp[3])+Cp*(vp[1]-vp[2])+(v[1]-v[2])/Rp+(v[1]-v[3])/Ru
	eqs[2] = Ce*vp[2]+Cp*(vp[2]-vp[1])+Cce*(vp[2]-vp[3])+(v[2]-v[1])/Rp+v[2]/Re+(v[2]-v[3])/Ro-gm*(v[1]-v[2])
	eqs[3] = Cu*(vp[3]-vp[1])+Cce*(vp[3]-vp[2])+Cc*(vp[3]-vp[4])+(v[3]-v[1])/Ru+(v[3]-v[2])/Ro+v[3]/Rc+gm*(v[1]-v[2])
	eqs[4] = Cc*(vp[4]-vp[3])+v[4]/Rl
	eqs[5] = Cs*(vp[5]-vp[1])+(v[5]-vs(t))/Rs
end

s1 = odeivs(f,{0,t2},{0,0,0,0,0})
plot(s1)
Vm,t1,t2 = -.1, 5.e-3, 10.e-3 -- Long pulse (5, 10 msec)
s2 = odeivs(f,{0,t2},{0,0,0,0,0})
plot(s2)
write_data("list10_28.dat",s1,s2)
