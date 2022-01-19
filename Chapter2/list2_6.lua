-- File list2.6 -- Examples of some functions in init.lua
-- Uncomment next line if needed
--require"init" -- Load basic Lua functions

-- List known objects and type (in _G table)
help() -- Same as whatis(_G)
-- List entries in math table -- names and types
whatis(math)
-- Probe object type of math.sin and math.pi
whatis(math.sin); whatis(math.pi)
-- Make sin a globally known function
sin = math.sin; whatis(sin)
-- Probe a string
whatis(_VERSION)
