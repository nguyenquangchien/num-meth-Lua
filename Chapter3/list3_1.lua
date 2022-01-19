--/* File list3_1.lua Successive substitution code */

x1 = 0 -- Initial guess
x2 = 0
for i=1,40 do 
	x1 = ((3-2*x1+4*x1^2)/5)^.25
	x2 = -((3-2*x2+4*x2^2)/5)^.25
	print(i,x1,i,x2)
end

