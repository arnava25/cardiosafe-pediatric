# FAERS cache provenance

Parquets built 8 June 2026 from all 40 FAERS quarterly ASCII releases,
2015q1 through 2024q4, downloaded from FDA.

demo.parquet   16,144,530 rows, all ages
drug.parquet
reac.parquet

Nine source zips were deleted after parquet construction to reclaim disk
(2018q1-q4, 2019q1, q2, q4, 2021q4, 2023q3). This does NOT affect the
parquets, which contain all 40 quarters. Verified by quarter-level row
counts.

To rebuild from scratch, re-download all 40 quarters from
https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html
