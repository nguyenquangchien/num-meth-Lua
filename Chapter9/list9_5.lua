-- /* File list9_5.lua */ Analysis of joint confidence limits

require"cltwodim"

c = {{},{},{}} -- C[] arrays
read_data('list9_3a.dat',c) -- Read first 3 columns of file
clim,xyb,cxx = cltwodim(c[1],c[3],90) -- 90% confidence limits
clim2,xyb2,cxx2 = cltwodim(c[1],c[3],95) -- 95% confidence limits
for i=1,2 do
	for j=1,2 do 
		print('c['..i..'], 95% limits = ',clim2[i][j],
		'90% limits = ',clim[i][j])
	end
end
scatterplot(c[1],c[3],xyb,cxx,xyb2,cxx2)
write_data('list9_5.dat',c,xyb,cxx,xyb2,cxx2) -- Save for plotting

