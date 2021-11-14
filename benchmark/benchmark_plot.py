
import matplotlib.pyplot as plt
import pandas as pd
import json
import numpy as np

with open("benchmark/benchmark_commits.json", "r") as f:
    dic_commits = json.load(f)

def plot_test(test_name, div_factor=1, units="ps"):
    dic_slice = {}
    dic_commits_keys = list(dic_commits.keys())
    for commit_hash in dic_commits_keys[-10:]:
        d = dic_commits[commit_hash]
        commit_date = d["commit_date"]
        label = commit_hash+"\n"+commit_date
        dic_slice[label] = np.array(d[test_name]) / div_factor
    df_slice = pd.DataFrame(dic_slice)
    plt.figure(figsize=(10,5))
    df_slice.boxplot(column=list(dic_slice.keys()))
    plt.ylabel(f"time ({units})")
    plt.savefig(f"benchmark\img\{test_name}.png")

plot_test("agg_dimer_small", div_factor=1e6, units="ms")
plot_test("agg_dimer_medium", div_factor=1e9, units="s")
plot_test("agg_dimer_big", div_factor=1e9, units="s")