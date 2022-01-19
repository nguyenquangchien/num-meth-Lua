-- /* File list9_19.lua */ -- Data fitting with multiple independent parameters

require"nlstsq"

-- Section 1. Input data file and define fitting function
infile = 'Nelson'
yd,x1d,x2d,yy,ycalc = {},{},{},{},{}
read_data(infile..'.txt',yd,x1d,x2d)
nd = #yd; for i=1,nd do yy[i] = math.log(yd[i]) end

ft = function(x,c) -- Define function to fit data
	return c[1] - c[2]*x[2]*math.exp(-c[3]*x[3]) - x[1]
end

c = {2,5e-9,-5e-2}; nc = #c
step={1.2,1.2,1.2}
--- Section 2. Perform data fit and print fitted parameters
del,err,nn = nlstsq({yy,x1d,x2d},fw,ft,c,actv,step) -- Call fiting, print max iterations 
print('Newton iterations, err =',nn,err)
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
for i=1,nd do ycalc[i] = newtonfc(ft,{yy[i],x1d[i],x2d[i]},c) end
write_data(infile..'.dat',x1d,x2d,yd,ycalc,yy)