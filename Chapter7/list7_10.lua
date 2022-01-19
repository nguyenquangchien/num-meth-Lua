-- /* File list7_10.lua */ -- First example of data fitting 
  
require"nlstsq"
 
xd,yd,i = {},{},1
for t=0,5,.1 do -- Generate some data
 	xd[i],yd[i] = t, 2*(1-math.exp(-.8*t))
  	i = i+1
end
  
ft = function(x,c) -- Define function to fit data
  	return c[1]*(1 - c[2]*math.exp(c[3]*x[2])) -x[1]
end
  
c = {4,.5,-1} -- Initial guess at coefficients
actv,step = {1,1,1},{0,0,0} -- Set up arrays
nd,nx,nc = #xd,2,3
fw = {} -- Weighting factors 
for i=1,nd do fw[i] = 1 end
x = {yd,xd} -- Data values
del,err,nmax = nlstsq(x,fw,ft,c,actv,step) -- Call fiting 
print(err,nmax) -- print results
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
