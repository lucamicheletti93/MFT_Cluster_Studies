# MFT_Cluster_Studies
Command to run MFTAssesment on data (ctf)
```ruby
o2-ctf-reader-workflow --onlyDet MFT --ctf-input o2_ctf_run00523399_orbit0776146286_tf0000616039_epn037.root | o2-mft-reco-workflow --shm-segment-size 15000000000 --clusters-from-upstream --disable-mc | o2-mft-assessment-workflow --disable-mc -b --run
```
The example works with:
/alice/data/2022/LHC22m/523399/raw/1230/o2_ctf_run00523399_orbit0776146286_tf0000616039_epn037.root
