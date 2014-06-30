from export import export_od
from testenv import cluster

ods = range(91, 1640+1)

for freq in 30, 44, 70:
    cluster.run(export_od, ods, common=[freq])
