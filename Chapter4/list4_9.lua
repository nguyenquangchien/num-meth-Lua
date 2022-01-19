-- File list4_9.lua  -- Biasing of a BJT with newton()

require"nsolv"; require"newton"

Is,vt,B,Va = 1.e-16,.026,100,100 -- Basic transistor parameters
Vcc,R1,R2,Rs,Rl = 15,20.e3,1.e3,200,2.e3 -- Basic circuit parameters
VC = 7.5 -- Desired collector voltage

fbjt = function(i,v)
	local ib,ic = v[4],v[5] -- Keep definitions local to this function
	i[1] = v[1]/R2 + (v[1] - Vcc)/R1 + ib
	i[2] = v[2]/Rs -  ib - ic
	i[3] = (v[3]-Vcc)/Rl + v[5]
	i[4] = ib - Is*(math.exp((v[1] - v[2])/vt)-1)
	i[5] = ic - B*ib*(1 + (v[3] - v[2])/Va)
end
v = {0,0,0,0,0} -- Initial guesses, try other values?
fc = function(R) -- Function for newton() call
	R2 = R; nsolv(fbjt,v) -- Use last values for v[]
	return v[3] - VC
end
R2 = newton(fc,R2)
printf('Calculated resistor value = %d\n',R2)
table.foreach(v,print)
