import DelayDiffEq # pkg_names_to_import
import DiffEqCallbacks # pkg_names_to_import
import LinearAlgebra # pkg_names_to_import
import OrdinaryDiffEq # pkg_names_to_import
import QuadGK # pkg_names_to_import
import QuantumOpticsBase # pkg_names_to_import
import Random # pkg_names_to_import
import Random # pkg_names_to_import
import Reexport # pkg_names_to_import
import ResumableFunctions # pkg_names_to_import
import SparseArrays # pkg_names_to_import
import StableRNGs # pkg_names_to_import

import Pkg
for (uuid, info) in Pkg.dependencies()
if info.name in ["DelayDiffEq", "DiffEqCallbacks", "LinearAlgebra", "OrdinaryDiffEq", "QuadGK", "QuantumOpticsBase", "Random", "Random", "Reexport", "ResumableFunctions", "SparseArrays", "StableRNGs"] # pkg_names_to_test
ENV["PREDICTMD_TEST_PLOTS"] = "true"
ENV["PREDICTMD_TEST_GROUP"] = "all"

include(joinpath(info.source, "test", "runtests.jl"))
end
end

