local intervalFile = arg[1]
local window = arg[2] or 0

for l in io.lines (intervalFile) do 
	local c, i = l:match("(%w+):(%d+)")
	print(c..":".. (i-window).. "-" .. (i+window))
end
