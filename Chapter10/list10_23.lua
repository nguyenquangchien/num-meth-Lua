--- /* File list10_22.lua */
-- Programs to integrate differential equation for pendulum
require"odeiv"
odebiv = odeb12

g,L,c,m = 983.21, 100, 400, 20
dc,w2 = c/(m*L), g/L
th0 = 0.95*math.pi

f1 = function(eqs,t,x,xp,xpp)
	eqs[1] = xpp[1] + dc*xp[1] + w2*math.sin(x[1])
	eqs[2] = xpp[2] + dc*xp[2] + w2*x[2] 
	eqs[3] = xpp[3] + w2*math.sin(x[3]) -- No Damping
	eqs[4] = xpp[4] + w2*x[4] -- No Damping
end

s1,err1 = odeivse(f1,{0,20},{th0,th0,th0,th0},{0,0,0,0})
plot(s1)
write_data("list10_23.dat",s1)