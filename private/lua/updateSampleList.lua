--[[
-- Updates a list of BAM files to the latest version in the picard aggregation path
-- Usage:
--
-- lua updateSampleList.lua samples.list > updated_samples.list
 ]]
function latestVersion(sample)
  local version = tonumber(sample:match("/v(%d+)/"))
  io.stderr:write(version .. " => ")
  f = io.open(sample)
  while (f == nil) do
    version = version + 1
    sample = sample:gsub("/v(%d+)/", "/v"..version.."/")
    f = io.open(sample)
  end
  io.stderr:write(version .. "\n")
  return(sample)
end

for sample in io.lines(arg[1]) do
  io.stderr:write("Updating sample: " .. sample .. " version: ")
  print(latestVersion(sample))
end
