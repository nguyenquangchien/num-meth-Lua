-- /* File find.lua */
-- Lua file for finding a string in files

local pat = {'%%','%(','%)','%[','%$','%?','%.','%+','%-','%*','%^'}
local ppat = {'%%%%','%%%(','%%%)','%%%[','%%%$','%%%?',
	'%%%.','%%%+','%%%-','%%%*','%%%^'}
string.nomagic = function(str)
	for i=1,11 do str = string.gsub(str,pat[i],ppat[i]) end
	return str
end
find = function(str)
	--str = string.nomagic(str)
	os.execute( "dir /-W /B "..filetype.." >$$tmpdir.txt")
	local ln
	for dname in io.lines('$$tmpdir.txt') do
		ln = 1
		for line in io.lines(dname) do
			--if string.find(line,str) then -- use below without nomagic()
			if string.find(line,str,1,true) then
				print(ln,dname,'\t',line)
			end
			ln = ln+1
		end
	end
	os.remove("$$tmpdir.txt")
end
setfenv(find,{filetype='*.lua',io=io,print=print,os=os,string=string})

getfenv(find).filetype = '*.lua'

find"nlstsq"
