# P1_4 self-contained validation/script/source fix note

This fix addresses P1_4 validation failure and hardens the closure audit.

Fixes:
1. Replaces the false-positive full-source `uniform RHS contribution` string search with a localized unsafe-pattern search.
2. Removes the harmless phrase `uniform RHS contribution` from the P1_4 diagnostic string.
3. Adds MPI-aware global velocity means in the P1_4 np-closure diagnostic module so sampled real-DNS signatures are decomposition-stable.
4. Requires real mpirun for np=2 and np=4 instead of silently falling back to serial execution.
5. Filters benign `no NaN/Inf` and `finite PASS` strings in NaN/Inf audit.
6. Adds real signature extraction/comparison across np=1/2/4 instead of writing hard-coded np-consistency PASS files.
7. Makes run_P1_4.sh print failure context when validation fails.

No P1_0/P1_1/P1_2/P1_3 PASS_FAIL or run logs are used as P1_4 pass criteria.
