-- File list2.1 -- Simple Lua operations

-- Test for maximum integer value 
function tst(n) -- function definition
	x = n or 1
	for i=1,400 do -- Do loop
		x = x + x -- double x value
		local y = x+1 -- Local variable
		z = i -- new variable
		if y-1~=x then -- Simple if statement
			return i -- Return from function
		end -- End if test
	end -- End for loop
end -- end function

i = tst() -- Call function with return value

 print("i, x, y, z = ",i,x,y,z) -- Simple print
 
-- Formatted print using io and string libraries
io.write(string.format('x = %40.10f \nx-1 =%40.10f \n',x,x-1))

-- Use of math library
x,y = math.sin(math.pi/4), math.cos(math.pi/4)
print('x, y = ',x,y)

-- Illustration of table
t = {n = 4, 1,3,5}
table.foreach(t,print) -- print table using library




