-- /* File listing6.9.lua */ -- Interpolation for several tables
require"intp"
require"spline"

--xd = {-2,-1.5,-1,-.5,0,.5,1,1.5,2} -- Example of pulse-like function
--yd = {0,0,0,.87,1,.87,0,0,0} -- example a
--xmin,xmax,dx = -2,2,.01
--xd = {1,2,3,4,5,6,7,8,9,10} -- Example of linear line with step
--yd = {3.5,3,2.5,2,1.5,-2.4,-2.8,-3.2,-3.6,-4}  -- example b
--xmin,xmax,dx = 1,10,.02
--xd = {-4,-3,-2,-1,0,1,2,3,4} -- Impulse like function
--yd = {0,0,0,0,1,0,0,0,0} -- example c
--xmin,xmax,dx = -4,4,.02
--xd = {-1,-.5,0,.5,1} --Second pulse like
--yd = {0.0385,0.1379,1,0.1379,0.0385} -- example d
--xmin,xmax,dx = -1,1,.005
--xd = {0,1,2,3} -- Sparse data
--yd = {0,1,4,3} -- example e
--xmin,xmax,dx = 0,3,.05
--xd = {40,48,56,64,72} -- vapor pressure of water
--yd = {55.3,83.7,123.8,179.2,254.5} -- example f
--xmin,xmax,dx = 40,72,.05
xd = {-200,-100,0,100,200,300,400,500} -- Type T thermocouple
yd = {-4.111,-2.559,-.67,1.517,3.967,6.647,9.525,12.575} -- g
xmin,xmax,dx = -200,500,1

flcb = intpf(xd,yd) -- Define local cubic interpolation fcn
fcsp = splinef(xd,yd) -- Define cubic spline interpolation fcn

x,ylcb,ycsp = {},{},{}
xv,i=xmin,1-- Calculate interpolated values
for xv=xmin,xmax+.001,dx do
	x[i], ylcb[i], ycsp[i] = xv, flcb(xv), fcsp(xv)
	i = i+1
end
write_data("list6_9table.dat",xd,yd)
write_data("list6_9interp.dat",x,ylcb,ycsp)






