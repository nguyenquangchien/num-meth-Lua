-- /* File list5_3.lua */

require"deriv"

f1 = function(x) return x^5 end
f2 = math.exp
val1 = 5*(1.e8)^4
val2 = math.exp(20)

val,n,ex = deriv(f1,1.e8)
print(val,n,ex,(val-val1)/val1)

val,n,ex = deriv(f2,20)
print(val,n,ex,(val-val2)/val2)

val,n,ex = rderiv(f2,20)
print(val,n,ex,(val-val2)/val2)

val,n,ex = lderiv(f2,20)
print(val,n,ex,(val-val2)/val2)

val,n,ex = deriv2(f1,1.e8)
val2 = 20*(1.e8)^3
print(val,n,ex,(val-val2)/val2)

f3 = math.sin
angle,vald,err = {},{},{}
for i=1,401 do
	ang = 2*math.pi*(i-1)/400
	val = math.cos(ang)
	angle[i] = ang
	vald[i] = deriv(f3,ang)
	if vald[i]==0 then err[i] = 0
	else err[i] = ((vald[i]-val)/val) end
end
write_data("list5_3.dat",angle,vald,err)

