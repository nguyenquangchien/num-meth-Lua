-- /* File list10_16.lua */
-- Programs to integrate first order diff. equation using odebrk()
require"odeiv"
odebiv = odebrk

f1 = function(eqs,t,y,yp) --  equation for cos(t), sin(t)
	eqs[1] = yp[1] - y[2]
	eqs[2] = yp[2] + y[1]
end
-- RK solution with one step interval, also error evaluation
s1,err1,nx,sx = odeive(f1,{0,100,2000},{1,0})
print(errstat(err1[2]))
s2,err2 = odeive(f1,{0,100,1000},{1,0})
print(errstat(err2[2]))
plot(s1[1],s1[2],err1[2])

write_data("list10_16.dat",s1,err1,s2,err2,sx)
