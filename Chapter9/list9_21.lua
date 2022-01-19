-- /* File list9_21.lua */ -- MOS transistor, improved model

require"nlstsq"

-- Section 1. Input data file and define fitting function
infile = 'mosiv9.20' 
id,vd,vg,vgd,ycalc = {},{},{},{},{}
read_data(infile..'.txt',id,vd,vg); nd = #id 

ft = function(x,c) -- Define function to fit data
	local vgt,vd = x[3] - c[2], x[2]
	local vdx = 0.5*(vd+vgt - math.sqrt((vd-vgt)^2))
	if vgt<0 then return -x[1] end
	return c[1]*(vgt*vdx - 0.5*vdx^2)*(1 + c[3]*vd)/(1+vgt*c[4]) - x[1]
end

c = {1.e-4,.20,0.01,0}; nc = #c
step = {0,0,5,0}

--- Section 2. Perform data fit and print fitted parameters
del,err,nn = nlstsq({id,vd,vg},fw,ft,c,actv,step) -- Call fiting, print max iterations 
for i=1,nc do printf('c[%d] = %12.4e +/- %12.4e\n',i,c[i],del[i]) end
for i=1,nd do ycalc[i] = newtonfc(ft,{id[i],vd[i],vg[i]},c) end
write_data(infile..'c.dat',id,vd,vg,ycalc)
plot(vd,id,ycalc)
