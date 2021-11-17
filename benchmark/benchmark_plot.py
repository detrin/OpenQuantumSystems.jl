import matplotlib.pyplot as plt
import pandas as pd
import json
import numpy as np

with open("benchmark/benchmark_commits.json", "r") as f:
    dic_commits = json.load(f)


def plot_test(test_name, div_factor=1, units="ps"):
    dic_slice = []
    dic_commits_keys = list(dic_commits.keys())
    plt.figure(figsize=(10, 5))
    plt.grid(True)
    # labels = dic_commits_keys[-10:]
    labels = []
    for commit_hash in dic_commits_keys[-10:]:
        d = dic_commits[commit_hash]
        if test_name in d:
            commit_date = d["commit_date"]
            label = commit_hash + "\n" + commit_date
            print(commit_hash, test_name, len(d[test_name]))
            dic_slice.append([x / div_factor for x in d[test_name]])
            # df_slice = pd.DataFrame(dic_slice)
            labels.append(label)
    plt.boxplot(dic_slice, showfliers=False)
    plt.xticks(range(1, len(labels) + 1), labels)
    plt.ylabel(f"time ({units})")
    plt.savefig(f"benchmark\img\{test_name}.png")


plot_test("agg_dimer_small", div_factor=1e6, units="ms")
plot_test("agg_dimer_medium", div_factor=1e6, units="ms")
plot_test("agg_building", div_factor=1e6, units="ms")
plot_test("getFranckCondonFactors", div_factor=1e6, units="ms")
plot_test("getFCProd", div_factor=1e9, units="s")
plot_test("getAggHamiltonian", div_factor=1e6, units="ms")
plot_test("getAggHamiltonianInteraction", div_factor=1e6, units="ms")

plot_test("trace_bath", div_factor=1e6, units="ms")
plot_test("evolutionExact", div_factor=1e6, units="ms")
plot_test("evolutionApproximate", div_factor=1e6, units="ms")
plot_test("schroedinger", div_factor=1e6, units="ms")
plot_test("liouvilleVonNeumann", div_factor=1e6, units="ms")
plot_test("evolutionOperatorIterator", div_factor=1e6, units="ms")
plot_test("master_int", div_factor=1e6, units="ms")
plot_test("master_ansatz", div_factor=1e6, units="ms")
