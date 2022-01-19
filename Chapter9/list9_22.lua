-- /* File list9_22.lua */ -- Response surface analysis

requires("nlstsq","prob")

-- Section 1. Input data file and define fitting function
infile = 'rspsurf' 
p,gr,y,ycalc,ydata = {},{},{},{},{}
read_data(infile..'.txt',{},p,gr,y) -- For uniformity analysis
--read_data(infile..'.txt',{},p,gr,{},y) -- For stress analysis
nd = #y

ft = function(x,c) -- Define function to fit data
	local x1,x2 = x[2],x[3]
	return c[1]+c[2]*x1+c[3]*x2+c[4]*x1^2+c[5]*x2^2+c[6]*x1*x2 - x[1]
end

c = {0,0,0,0,0,0}; nc = #c
--actv = {1,1,1,0,0,1} -- Uncomment for eliminating c[4] and c[5]

--- Section 2. Perform data fit and print fitted parameters
del,err,nn = nlstsq({y,p,gr},fw,ft,c,actv) 
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e \t (+/-ratio)%9.4f\n',
	i,c[i],del[i],del[i]/math.abs(c[i])) end
x1,x2,i = {},{},1
for xx1=4,80.1,1 do -- Points for surface plot
	for xx2=2,9,.1 do
		x1[i],x2[i],ycalc[i],i = xx1,xx2,ft({0,xx1,xx2},c),i+1
	end
end
for i=1,nd do ydata[i] = newtonfc(ft,{y[i],p[i],gr[i]},c) end
_,_,rsq = clinear(ydata,y)
print('Rsquared for model = ',rsq)
write_data(infile..'a.dat',y,p,gr,ydata,{ycalc,x1,x2})

