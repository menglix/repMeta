# repMeta (development version 0.1.1)

## Output naming changes

* Renamed the `out_studies` element returned by `Rm.func.iterative()` and
  `Rm.func.iterative.boot()` to `nonrep_studies`, to match the terminology
  ("non-replicable studies") already used everywhere else in the package's
  documentation and examples. Code doing `result$out_studies` needs to be
  updated to `result$nonrep_studies`.

## Bug fixes

* Fixed `to.dat.repMeta()` throwing an error under `metafor` >= 3.8.1:

  ```
  "Error in escalc(measure, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, :
Unknown 'to' argument specified.",
  ```

  **Cause:** `repMeta` was built against `metafor` 2.4.0, where `escalc()` accepted values passed positionally. Later updates in `metafor` redesigned `escalc()` to resolve its own arguments via non-standard evaluation, and now  generates errors whenever those optional arguments were left unspecified.

  **Fix:** `to.dat.repMeta()` now builds the `escalc()` call as a named argument list, drops any `NULL` entries before dispatching, and handles the `digits` argument via `missing()` so it's never force-evaluated when unsupplied.

The fix produces byte-identical output to the original implementation when run
against the `metafor` version `repMeta` was originally built for, and resolves
the failure on current `metafor`.

# repMeta 0.1.0

* Initial release at 03June2023.
