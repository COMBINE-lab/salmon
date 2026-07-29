# salmon 2.4.1

A bug-fix release on 2.4.0. It fixes a race in the RAD reader that could cause
`salmon quant --rad` and alignment mode to **silently produce zero or partial
counts** — an empty or short `quant.sf`, with no error and no warning.

No index rebuild is required, and no other behavior changes. If your 2.4.0 runs
produced sensible output, they were unaffected; the failure mode is conspicuous
rather than subtle when it occurs.

## Fix: RAD reader workers could exit with chunks still queued (#1071)

`libradicl`'s reader pushes every meta-chunk onto a shared queue and only then
sets its `done` flag. Salmon's worker threads exited on:

```rust
if let Some(mc) = queue.pop() { … }
else if done.load(Ordering::Acquire) { break }
```

A worker whose `pop()` came up empty in the window *just before* the producer
pushed would then observe `done` and break, abandoning everything still queued.
Four worker loops had this shape. Each now takes one more pass over the queue
after first observing `done`, which closes the window.

### How often could this happen?

Rarely, but silently, which is why it is worth patching promptly:

- Observed at roughly **1 in 600** runs under deliberate CI stress on aarch64
  (the whole test binary running concurrently, which is what made worker threads
  slow enough to lose the race).
- **Never observed on x86_64** outside artificially forced conditions.
- Most likely at **low thread counts**. With a single worker there is no sibling
  to drain the queue, so one lost race loses the entire file. With the window
  held open artificially, however, every worker can lose it simultaneously
  regardless of thread count — so a higher `-p` reduces the odds rather than
  eliminating them.

Affected paths are `salmon quant --rad`, alignment mode (`-a`), and the
`--deterministic` intermediate requant, all of which read RAD/BAM through the
parallel reader. The reads-based mapping path does not use it.

If you have 2.4.0 results you are unsure about, the signature is unambiguous:
`num_processed` in `aux_info/meta_info.json` will be short of the fragments the
input actually contains, and in the worst case `quant.sf` is all zeros.

## A note on how this was found, and on the test that missed it

The bug surfaced as an intermittent CI failure of the `thread_count_invariant`
test, which compares three quantification runs at different thread counts. That
framing was misleading in two ways worth recording.

It is **not a thread-count-dependent bug**. The instrumented reproduction showed
the 1-thread and 8-thread runs agreeing exactly, and the *repeat* of the 1-thread
run returning zero. Any worker can lose the race.

More importantly, the test **could not have caught the worst case**. It asserts
only that the three runs agree with each other, so a fault that affects every run
looks like a pass — with the race window forced open, all three runs returned
zero and the test reported success. The intermittent CI failures were the *lucky*
case, where only some runs lost the race.

The regression test added with the fix therefore asserts absolute correctness —
that every fragment in the input is processed and mass is conserved — rather than
mutual consistency, and it was verified to fail against the unfixed code.

## Upstream

`libradicl` itself is not affected: it only ever produces into the queue and
never consumes, so the library contains no instance of this. Its two parallel
examples do use the racy shape, and are being reported upstream along with a
proposal for a drain-safe iterator, so that consumers do not each have to
rediscover the ordering contract. `alevin-fry` guards correctly at all of its
call sites.

## Upgrading

Drop-in from 2.4.0: no index rebuild, no option changes, and results are
unchanged for any run that was not hit by the race.
