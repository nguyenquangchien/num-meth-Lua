-- /* File list8_13.lua */

require"prob"

redwell,whitney = {},{} -- Two data sets for comparison
read_data('redwell.txt',redwell)
read_data('whitney.txt',whitney)
plot({makeCDF(whitney)},{makeCDF(redwell)})

KS_test(redwell,whitney) -- K-S Distribution Test
