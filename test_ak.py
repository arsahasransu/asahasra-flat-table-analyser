import awkward as ak

arr = ak.Array([[1, 2, 3], [], [4, 5]])
pt = ak.Array([[10, 20, 30], [], [40, 50]])

mask = pt > 15
filtered = arr[mask]
filtered_non_empty = filtered[ak.num(filtered) > 0]

print(filtered)
print(filtered_non_empty)
