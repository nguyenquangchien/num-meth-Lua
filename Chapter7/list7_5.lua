--/* File list7_5.lua */ Fourier analysis of sunspot activity

require"Fourier"

jk,yyd,ssd = {},{},{}
read_data('monthssn.dat',jk,yyd,ssd) -- read data
nmths = #yyd -- Number of months
plot(yyd,ssd) -- Look at data
fx,fy,fz = fourier(ssd,100) -- Get Fourier components
freq,fmag = {},{} -- Tables of magnitude and angle
for i=1,#fx do 
	freq[i] = nmths/(12*i) -- Convert from harmonic nnbr to frequency
	fmag[i] = fx[i][1]
end

plotFourier(fx) -- View Fourier components
plot(iFourier(fx)); plot(freq,fmag)
write_data('list7_5.dat',freq,fmag) -- Save analysis