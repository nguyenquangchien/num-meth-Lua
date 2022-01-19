-- /* File elemfunc_test.lua */ -- Some tests of elementary functions
-- Probability function, Error function, Complementary Error function, Gamma function

require"elemfunc"

--print(elemfunc.E(1-.9e-14),elemfunc.E(1-1.1e-14))
i,x2t,y2t,y3t,dt,y4t = 1,{},{},{},{},{} -- Now test function
local c=os.clock()
i = 1
for x=-5,5.0001,.01 do 
--for x=0.0,4.0001,.1 do
if x~=0 then 

	x2t[i] = x
	y2t[i] = elemfunc.gamma(x)
	--y2t[i] = prob(x)
	--y2t[i] = en(x,1)
	--y2t[i] = elemfunc.Ci(x)
	y3t[i] = elemfunc.Ci(x)
	y4t[i] = elemfunc.Cin(x)
	--y2t[i] = K(x)
	--y3t[i] = E(x)
	--y2t[i] = elemfunc.Pn(x)
	--y3t[i] = elemfunc.Qn(x)
	--y4t[i] = elemfunc.erf(x)
	--y2t[i] = elemfunc.Jn(x,3)
	--y3t[i] = elemfunc.Yn(x,3)
	i = i+1
end	
end
print('time =',os.clock()-c)
plot(x2t,y3t) ;plot(x2t,y4t)
--plot(x2t,y2t); plot(x2t,y3t)
write_data("test.dat",x2t,y2t,y3t,y4t)
--write_data("test.dat",x2t,y2t,y3t)
