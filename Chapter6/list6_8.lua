-- /* File list6_8.lua */ -- Interpolation example with amplifier data
require"intp"

xd,yd,_ = {},{},{} -- dummy table _ = {}
read_data('list4_10.dat',yd,_,_,_,xd)

npt,xm,im = #xd, 0, 0 
-- Search for largest value in gain -- xd here
for i=1,npt do if xd[i]>xm then xm,im = xd[i],i end end

h3db = xm-10*math.log10(2) -- -3Db gain value
fl = intp(xd,yd,h3db) -- Lower -3dB frequency

i=1
for j= npt,im,-1 do -- Define table for upper -3dB range
	 xd[i],yd[i],i = xd[j], yd[j], i+1
 end -- Reverse x data for increasing table values
 
fh = intp(xd,yd,h3db) -- Upper -3dB frequency
printf("Maxumim gain in dB =   %12.4e dB\n",xm)
printf("Lower 3dB Frequency = %12.4e Hz\n",fl)
printf("Upper 3dB Frequenty = %12.4e Hz\n",fh)
printf("Amplifier Bandwidth   = %12.4e Hz\n",fh-fl)
