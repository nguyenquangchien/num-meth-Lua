--- /* File list10_27.lua */
-- Programs to integrate surface charge differential equation

require"odeiv"

a,a0 = {.1,.2,.5,1,2,5,10}, {0,0,0,0,0,0,0}
na = #a

f = function(eqs,t,f,fp) -- Define equations with multiple a values
	for i=1,na do eqs[i] = f[i]*fp[i] + f[i]*(1-a[i]) - a[i] end
end

s1,err1 = odeivse(f,{0,100},a0)
print(unpack(maxvalue(s1))) ; print(unpack(maxvalue(err1)))
plot(s1); plot(err1)

write_data("list10_27.dat",s1,err1)

