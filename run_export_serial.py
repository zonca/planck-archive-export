from export import export_od
from testenv import cluster

ods = range(91, 1640+1)
ods = range(219, 269+1)

freqs = 44, 70
freqs = 70,

for freq in freqs:
    for od in ods:
        export_od(od, freq)
