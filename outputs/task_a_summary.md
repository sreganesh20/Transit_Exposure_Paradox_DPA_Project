# Task A: Geographic Coverage Validation Summary

- **Total CTA square queries evaluated:** 426
- **Fully covered queries:** 340 (79.81%)
- **Partially covered queries:** 36 (8.45%)
- **Outside Chicago queries:** 50 (11.74%)

## Zero-Crime Square Queries (`crimes_2024 == 0`)

- **Total zero-crime square queries:** 54
- **Outside Chicago:** 50
- **Partially covered:** 4
- **Fully covered:** 0

## Key Takeaway

All 54 zero-crime square queries have incomplete geographic coverage (50 outside Chicago and 4 partially covered), indicating that jurisdictional boundaries likely contribute substantially to the apparent prediction errors. For queries completely outside Chicago, the absence of Chicago police records is expected; for partially covered queries, incomplete overlap indicates boundary exposure rather than proving clipping alone caused a zero count.
