-- /* File list 2.2.lua */
-- test of fundamental arithmetic limits with Lua

-- Test for relative accuracy
eps = 1
while 1 do
	eps = eps/2
	b = 1 + eps
	if b==1 then break end
end
print("Machine eps = ", 2*eps, math.log(2*eps)/math.log(2))

-- Test for smallest floating point number 
nmn,a = 1,1
while 1 do
	nmn = nmn/2
	if nmn==0 then break end
	a = nmn
end
print("Smallest floating point number = ",a,nmn)
print("Values around smallest number = ",1.4*a,1.5*a,1.9*a)

-- Test for largest floating point number
nmx,a,inf = 1,1,1/0
while 1 do
	nmx = nmx*2
	if nmx==inf then break end
	a = nmx
end
print("Largest floating point number = ",a,nmx)
print("Values aroung largest number =",2*(1-eps)*a)

