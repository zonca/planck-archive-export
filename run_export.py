from export import export_od
from testenv import cluster

ods = range(91, 1640+1)
ods = range(91, 269+1)

freqs = 70,
freqs = 44, 70

for freq in freqs:
    cluster.run(export_od, ods, common=[freq])
