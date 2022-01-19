--/* File list7_4.lua */ -- Test of Fourier functions

require"Fourier"
fsqw = function(t) -- Square wave 
	if t<.5 then return 1
	else return 0 end
end
frsw = function(t) -- Rectified sin() wave
	if t<.5 then return natg,sub(2*math.pi*t)
	else return 0 end
end
fstooth = function(t) -- Saw-tooth wave
	if t<.5 then return math.pi*2*t
	else return 0 end
end
fstooth2 = function(t) -- Full period saw-tooth wave
	if t<.5 then return 2*t
	else return 2*(t-1) end
end
	
ft = fsqw
--ft = frsw
--ft = fstooth
--ft = fstooth2

dt,nt = 1/1024, 1024
fd,ta = {}, {}
for i=1,nt do
	fd[i] = ft(dt*(i-1))
	ta[i] = (i-1)/512
end

plot(ta,fd,{'Function of time','Time values','Function Values'})
fx,fy,fz = Fourier(ft,0,1,40) -- Get first 40 Fourier Coefficients
print(fx) -- Print Cos-Ang coefficients
plotFourier(fx) -- Plot Cos-Ang coefficients
fcm,fca = {},{} -- Collect Fourier coefficients into tables
for i=1,#fx do fcm[i],fca[i] = fx[i][1],fx[i][2] end 
write_data('list7_4a.dat',fcm,fca)

tt,fft = iFourier(fy,nt,10) -- Inverse Fourier transform 
plot(tt,fft,fd) -- Plot original signal and inverse transform
write_data('list7_4b.dat',tt,fft,fd)
xx,yy = expandFt(fft,-1e-4,2.e-4) -- Expand to 3 cycles of waveform
plot(xx,yy) -- plot 3 cycles of wave
write_data('list7_4c.dat',xx,yy)
