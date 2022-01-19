-- /* File list8_14.lua */ MC analysis of mean and std limits

require"prob"
local Pn = elemfunc.Pn

yield = {}
read_data('steel_yield.txt',yield) -- Read data file -- change as desired
nd = #yield; snd = math.sqrt(nd)
mean,std = stats(yield) -- Mean and std
xd,yd = makeCDF(yield) -- CDF of sample data
print('Number of data, mean, std = ',nd,engr(mean),engr(std))
nmc = 10000 -- Number of Monte Carlo runs

Pnms = function(x,m,s) return Pn((x-m)/s) end -- Gaussian with m and s

ymn,xmn = {},{}
mmc,mstd,Dn = {},{},{} -- collect mean and std for MC runs
for j=1,nmc do -- Begin MC runs
	for i=1,nd do -- Generate MC data set
		ymn[i] = std*gnormal()+mean
	end
	mmc[j],mstd[j]	= stats(ymn) -- Now evaluate mean and std
	xmn,ymn = makeCDF(ymn) -- CDF for Dn evaluation
	Dn[j] = 0
	for i=1,nd do -- Step over data set, collect Dn data
		dt = math.abs(ymn[i] - Pnms(xmn[i],mean,std))
		if dt>Dn[j] then Dn[j] = dt end
	end
end
pbv = {68.27,95.45,90,95,99,99.9}; np = #pbv
for j=1,nmc do Dn[j] = snd*Dn[j] end -- Final K-S statistic
x,y = makeCDF(mmc) -- CDF of MC mean values
for i=1,np do print('Limits for mean at '..engr(pbv[i],'%')..' = ',climits(x,pbv[i])) end
write_data('list8_14a.dat',x,y)
x,y = makeCDF(mstd) -- CDF of MC std values
for i=1,np do print('Limits for std at     '..engr(pbv[i],'%')..' = ',climits(x,pbv[i])) end
write_data('list8_14b.dat',x,y)
x,y = makeCDF(Dn) -- CDF of MC D values
for j=1,nmc do y[j] = 1 - y[j] end
write_data('list8_14c.dat',x,y)

