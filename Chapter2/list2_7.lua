-- File list2.7 -- Examples of reading and writing tabluar data

x,y,n = {},{},1000
for i=1,n do -- Define tables of values
	x[i] = (i-1)/(n-1)
	y[i] = math.sin(6*math.pi*x[i])
end
plot(x,y) -- pop-up plot
write_data('list2_7.dat', x,y) -- Save data to disk file

xy = {{},{}} -- New table of tables for reading data
read_data('list2_7.dat',xy) -- Retreive two column data
plot(xy[1],xy[2]) -- pop-up plot of retreived data
-- xy[1] is same as x and xy[2] is same as y data
whatis(xy)
