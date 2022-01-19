-- File list2.1 -- Simple Lua operations

-- Test for maximum integer value 
sum = 1;x = .5;fact = 1
val = math.exp(.5);print('val =',val)
xx = {{},{},{},{},{},{}}
for i=1,20 do
	fact = fact*i
	trm = x^i/fact
	sum = sum+trm
	xx[1][i],xx[2][i],xx[3][i],xx[4][i],xx[5][i],xx[6][i] =
	i,fact,trm,sum,trm/sum,(val-sum)/val
	print(i,fact,trm,sum,trm/sum,(val-sum)/val)
end
write_data('fig2.1.dat',xx)



