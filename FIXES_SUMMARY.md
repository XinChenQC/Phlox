# Parallel Processing Fixes for Phlox

## Changes Made

### 1. Fixed LAMMPS `skip_frames()` Performance (trajectory_reader.py)
**Problem:** Workers were calling `read_frame()` to skip, which fully parses each frame (reads all lines, creates Atom objects, sorts atoms). For Worker 7 skipping 6000 frames, this meant ~5 minutes of wasted parsing.

**Solution:** Implemented fast line-based skipping like ReaxANA:
- Cache `lines_per_frame` on first read
- Skip by reading raw lines without parsing: `readline()` × (n_frames × lines_per_frame)
- **Expected speedup:** ~100x faster skipping (from ~300s to ~3s for 6000 frames)

### 2. Reduced Worker Stagger Delay (analyzer.py)
**Problem:** 30-second delay per worker made debugging impossible:
- Worker 7 waited 210 seconds before starting
- Made it appear that workers were running serially

**Solution:** Reduced stagger from 30s to 2s
- Still prevents I/O contention on file open
- Workers now start much closer together
- **Total delay reduction:** 210s → 14s for 8 workers

### 3. Added Skip Progress Logging (analyzer.py)
**Problem:** Workers were silent during skip phase, appearing "stuck"

**Solution:** Added INFO-level logging:
```
Worker 1 will skip 893 frames before processing
Worker 1 skipping 893 frames...
Worker 1 skip complete, starting processing
```

### 4. Added Multiprocessing Guard (test_real_data.py)
**Problem:** Without `if __name__ == '__main__':`, worker processes re-execute the script, potentially causing issues

**Solution:** Wrapped all execution code in proper guard
- Ensures clean multiprocessing on all platforms
- Changed default to `n_cores=8` for testing

## Expected Performance Improvement

**Before:**
- 8 cores slower than 1 core (due to 3.5x redundant I/O)
- Total time: ~710s for 7151 frames

**After:**
- Skip time reduced from ~1750s total to ~17s total
- Stagger delay reduced from 210s to 14s
- **Expected speedup:** 6-7x faster than before (actual parallel speedup)
- Should now be ~8x faster than 1 core

## Test It

Run with 8 cores:
```bash
uv run python test_real_data.py
```

You should now see:
- All workers starting within 14 seconds
- Skip progress messages
- Interleaved [FRAG] output showing true parallelism
