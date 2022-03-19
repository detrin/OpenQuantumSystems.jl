_simplecontainergenerator_containers_temp_depot_tmpdir = mktempdir()
atexit(() -> rm(_simplecontainergenerator_containers_temp_depot_tmpdir; force = true, recursive = true))
_simplecontainergenerator_containers_temp_depot = joinpath(_simplecontainergenerator_containers_temp_depot_tmpdir, "simplecontainergenerator_containers_temp_depot")
mkpath(_simplecontainergenerator_containers_temp_depot)
pushfirst!(Base.DEPOT_PATH, "/opt/simplecontainergenerator_containers/julia_depot")
if get(ENV, "SIMPLECONTAINERGENERATOR_CONTAINER_NO_TEMP_DEPOT", "") != "true"
pushfirst!(Base.DEPOT_PATH, _simplecontainergenerator_containers_temp_depot)
end

