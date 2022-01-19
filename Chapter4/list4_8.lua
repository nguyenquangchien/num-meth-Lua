-- File list4_8.lua  -- Program to analyze the biasing of a BJT 

require"nsolv"

Is,vt,B,Va = 1.e-16,.026,100,100 -- Basic transistor parameters
Vcc,R1,R2,Rs,Rl = 15,20.e3,1.e3,200,2.e3 -- Basic circuit parameters

fbjt = function(i,v)
	local ib,ic = v[4],v[5] -- Keep definitions local to this function
	i[1] = v[1]/R2 + (v[1] - Vcc)/R1 + ib
	i[2] = v[2]/Rs -  ib - ic
	i[3] = (v[3]-Vcc)/Rl + v[5]
	i[4] = ib - Is*(math.exp((v[1] - v[2])/vt)-1)
	i[5] = ic - B*ib*(1 + (v[3] - v[2])/Va)
end
v = {0,0,0,0,0} -- Initial guesses
Ra,vb,ve,vc,ib,ic = {},{},{},{},{},{} -- Tables to store results
R2 = .1e3 -- Initial resistor value, zero causes problems with eqs
i = 1
while R2<=4e3 do
	nit,err = nsolv(fbjt,v) -- Use previous solution as guess
	if nit==NXX then print('Convergence error at i = ',i,err) end
	print(R2,v[3],v[5])		
	vb[i],ve[i],vc[i] = v[1], v[2], v[3]
	ib[i],ic[i],Ra[i] = v[4], v[5], R2
	i = i+1; R2 = R2 + .1e3
end
write_data("list4_8.dat",Ra,vb,ve,vc,ib,ic)
plot(Ra,vc,ve)