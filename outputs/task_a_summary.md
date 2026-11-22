# Task A: Geographic Coverage Validation Summary

- **Total CTA square queries:** 426
- **Fully covered:** 340
- **Partially covered:** 36
- **Outside Chicago:** 50
- **Zero-crime queries (`reported_crimes_2024 == 0`):** 57
  - Outside Chicago: 50
  - Partially covered: 7
  - Fully covered: 0

## Interpretation
- Queries classified as **outside Chicago** or **partially covered** have missing or clipped Chicago crime coverage where query footprints extend beyond city limits.
- Zero-crime queries that are **fully covered** inside Chicago are not explained by boundary truncation and reflect other factors (such as low crime counts within the query window/radius or station-level differences).
