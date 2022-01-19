-- /* File Listing7_6.lua */ -- Signal in noise example

require"fourier"

tpn,tsn = {},{} -- Read signal_noise and signal, 256 points
read_data('target_plus_noise.txt',tpn)
read_data('target_signal.txt',tsn)
fx,fy,fz = fourier(tpn,128) -- Fourier analysize
fm,fang = {},{} -- Magnitude and angle of Fourier comp
for i=1,#fx do fm[i],fang[i] = fx[i][1],fx[i][2] end
plotFourier(fx) -- Plot coefficients
ttn = {} -- Generate time points 
for i=1,#tpn do ttn[i] = i-1 end
plot(ttn,tpn) -- Plot signal plus noise
plot(ttn,tsn) -- Plot signal
tt,sig = iFourier(fx,256,60) -- Recover signal
write_data('list7_6a.dat',fm,fang)
plot(tt,sig) -- Plot recovered signal
write_data('list7_6b.dat',tt,sig,tpn,tsn)