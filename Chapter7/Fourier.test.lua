require"Complex"
require"Fourier"
--
ftx = function(t)
	--return math.cos(2*math.pi*t) + 2*math.sin(6*math.pi*t+math.pi/4)
	--return math.cos(2*math.pi*t)^2
	if t<.5 then return math.sin(2*math.pi*t)
	else return 0 end
	--if t<.5 then return 4*t
	--else return -4*(t-.5) end
end

fsqw = function(t) -- Square wave 
	if t<.5 then return 1
	else return 0 end
end
fsqw2 = function(t) -- Second square wave
	if t<.25 then return 1
	elseif t<.75 then return -1
	else return 1 end
end
fsqw3 = function(t)
	if t<.5 then return 1
	else return -1 end
end
frsw = function(t) -- Rectified Cos wave
	if t<.25 then return math.cos(2*math.pi*t)
	elseif t<.75 then return 0
	else return math.sin(2*math.pi*(t-.75)) end
end
frsw1 = function(t)
	if t<.5 then return math.sin(2*math.pi*t)
	else return 0 end
end
frsw2 = function(t)
	if t<.4 then return math.sin(2*math.pi*(t+.1))
	else return 0 end
end

fstooth = function(t)
	if t<.5 then return math.pi*2*t
	else return 0 end
end

fstooth2 = function(t)
	if t<.5 then return 2*t
	else return 2*(t-1) end
end
fstooth3 = function(t)
	if t<1 then return math.pi*t
	else return math.pi*(t-2) end
end


f910 = function(t)
	if t<.25 then return 1
	elseif t>.75 then return 1
	else return 0 end
end
	
--ft = fsqw	
--ft = fsqw2
--ft = fsqw3
--ft = frsw
--ft = frsw1
--ft = frsw2
--ft = fstooth
ft = fstooth2
--ft = fstooth3
--
--dt = 1/1024
dt = 1/512
fd = {}
ta = {}
--for i=1,1024 do 
for i=1,512 do
	fd[i] = ft(dt*(i-1))
		--ta[i] = (i-1)/1024
	ta[i] = (i-1)/512
end
--
plot(ta,fd,{'Function of time','Time values','Function Values'})

t = os.time()
fx,fy,fz = Fourier(ft,0,1,20)
print(fx)
print(fy)
print(fz)
--fx,fy,fz = Fourier(ft,-.1,.9,40)
print(type(fx),type(fy),type(fz))
plotFourier(fx)
plotFourier(fy)
plotFourier(fz)
--fx,fy,fz = Fourier(ft,0,1,50)

t1 = os.time()
print('time = ',t1-t)

tt,fft = iFourier(fy,512,20) --,512,3) --,1024)
--for i=1,table.getn(tt) do print(i,tt[i],fft[i]) end

print('time2 = ',os.time()-t1)
plot(tt,fft,fd)
xx,yy = expandFt(fft,-1e-4,2.e-4)
plot(xx,yy)
--plotcossin(fz)
--plotcostha(fy)
--plotexpang(fx)
--[[

jk,yyd,ssd = {},{},{}
read_data('monthssn.dat',jk,yyd,ssd)
nmths = table.getn(yyd)
table.foreach(ssd,print)
plot(yyd,ssd)
fx,fy,fz = fourier(ssd,100)
freq,fmag = {},{}
for i=1,table.getn(fx) do 
	freq[i] = nmths/(12*i)
	fmag[i] = fx[i][1]
end

plotFourier(fx)
plot(iFourier(fx))
write_data('list7.4.dat',freq,fmag)]]