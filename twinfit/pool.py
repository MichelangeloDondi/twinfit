"""A function over a list of jobs, in order, in several forked processes or in this one.

The jobs are not sent to the processes: a model that was wrapped by a mechanism, or a law written as a function inside
another, cannot be pickled. Each forked process reads the list that this module held when it was forked, and is sent the
place of its job alone. So the platform must be able to fork a process: where it cannot, the jobs run in this one.

Each process should run its linear algebra on one thread, or the processes fight for the cores: set
VECLIB_MAXIMUM_THREADS, OPENBLAS_NUM_THREADS or OMP_NUM_THREADS to 1 before numpy is first loaded.

A forked process is not safe on every platform, and the call never waits for a worker that is gone. On one platform a
forked process ended at its first decomposition of a matrix of 80 rows or more, once the calling process had made one
of that size: its library of linear algebra hands such a matrix to several threads, and a forked process cannot use
them. The pool that this module used until then waited for the lost rows and never came back. Now each worker asks
that library for one thread before its first job, where the library has a call for it (`_one_thread`), and a pool
whose worker ended with no result is closed: the rows that are missing are computed in the calling process, one after
another, with a warning, and `LOST` keeps the event. The rows are the same and the time is longer. Where the platform
cannot fork at all the same is done. What this does not cover: a job that ends its own process whoever runs it ends
the calling process too.

A job that raises in a worker is not sent back as its exception: one whose class the calling process cannot rebuild
from its record broke the pool, and was told as a lost worker. The worker gives back the place of the job. The jobs
that no worker has taken are dropped, those in flight are waited for, and the calling process runs that job again, so
that its exception comes out as it is. A refusal of a job therefore costs that job twice and up to one more job for
each worker. If the job does not raise again in the calling process, the fault was the worker's alone: the rows that
are missing are computed in the calling process, with the warning and the entry in `LOST` of a lost worker. What this
does not cover: a job that is interrupted from the keyboard in a worker, which is left to the pool. The warning and
the entry in `LOST` are written before the calling process computes the missing rows, so that a long wait has its
reason on the screen: if one of those rows then raises, the entry counts rows that were never computed.

A twin that is too long for one sitting can be read elsewhere, one synthetic set at a time, and put together later by
the same function of the package, so that no program outside the package repeats its arithmetic:

    with reading(one(17)) as got:      the function of the package reads its 18th set alone, and is left there
        admission(problem, result, sets=200)
    row = got.row                      what that set gave, to be stored

    with reading(stored(rows)):        the same call, later: its sets are the stored rows, in order, and nothing is fitted
        table = admission(problem, result, sets=200)

The jobs of such a function are made when they are read (Lazy), so the 18th set is drawn without the 17 before it.
"""
from __future__ import annotations

import sys
import warnings
from contextlib import contextmanager

_JOBS = []           # what a forked process reads
_READER = None       # how the sets of a twin are read, when a program reads them elsewhere
LOST = []            # each pool whose workers were lost: why, and how many rows the calling process then computed


class Lazy:
    """A list whose entries are made when they are read: `make(i)` gives entry i."""
    def __init__(self, make, n: int):
        self._make, self._n = make, int(n)

    def __len__(self) -> int:
        return self._n

    def __iter__(self):
        return (self._make(i) for i in range(self._n))

    def __getitem__(self, i: int):
        if not 0 <= i < self._n:
            raise IndexError(f"no entry {i} in a list of {self._n}")
        return self._make(i)


class _Left(Exception):
    """Raised by the reader of one set once that set is read: the function of the package is left before it puts a table together."""


class _Got:
    row = None
    read = False


def elsewhere() -> bool:
    """Whether a reader is in place: a function that reads its sets one after another by itself must then pass through `in_order`."""
    return _READER is not None


@contextmanager
def reading(reader):
    """Inside the block, every `in_order` of the package is done by `reader`. One block at a time."""
    global _READER
    if _READER is not None:
        raise RuntimeError("a reader is already in place: one block at a time")
    got = _Got()
    _READER = lambda fn, jobs: reader(fn, jobs, got)       # noqa: E731
    try:
        yield got
    except _Left:
        pass
    finally:
        _READER = None


def one(index: int):
    """The reader that runs the job at one place of the list and leaves. The other jobs are neither made nor run."""
    def reader(fn, jobs, got):
        if not 0 <= index < len(jobs):
            raise IndexError(f"the list of {len(jobs)} jobs has no place {index}")
        got.row, got.read = fn(jobs[index]), True
        raise _Left
    return reader


def stored(rows: list):
    """The reader that gives back rows that were read before, in the order of the jobs. A list of another length is refused:
    a table is put together from every set that it was asked for."""
    def reader(fn, jobs, got):
        if len(jobs) != len(rows):
            raise ValueError(f"{len(rows)} stored rows for {len(jobs)} sets: a table is put together from every set that it was asked for")
        got.read = True
        return list(rows)
    return reader


class _Raised:
    """What a worker gives back for a job that raised: the place of the job, which the calling process runs again."""
    def __init__(self, index: int):
        self.index = int(index)


def _run(index: int):
    fn, args = _JOBS[index]
    try:
        return fn(args)
    except KeyboardInterrupt:
        raise
    except BaseException:                                  # noqa: BLE001  whatever it is, it is raised again where the caller can read it
        return _Raised(index)


def _one_thread() -> None:
    """Run by each forked process before its first job: one thread for the linear algebra, where its library has a call
    for it. A forked process holds the one thread that forked it, and a library that hands a large matrix to threads
    of its own then ends the process. Where the call is absent nothing is done, and a lost worker is met by `in_order`."""
    if sys.platform != "darwin":
        return
    try:
        import ctypes
        call = getattr(ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate"), "BLASSetThreading", None)
        if call is not None:
            call.argtypes, call.restype = [ctypes.c_int], ctypes.c_int
            call(1)                                        # 1 is the single threaded setting of that library
    except Exception:                                      # noqa: BLE001  a platform without that library: nothing to ask
        pass


def in_order(fn, jobs: list, workers: int = 1) -> list:
    """fn(job) for each job, in the order of the list, whatever the number of workers. A pool whose worker ended with
    no result is closed and the missing rows are computed here, with a warning: the call does not wait for them."""
    if _READER is not None:
        return _READER(fn, jobs if hasattr(jobs, "__getitem__") else list(jobs))
    jobs = list(jobs)
    if workers <= 1 or len(jobs) < 2:
        return [fn(j) for j in jobs]
    global _JOBS
    _JOBS = [(fn, j) for j in jobs]
    rows, why, raised = [], "", False
    try:
        try:
            import multiprocessing
            from concurrent.futures import ProcessPoolExecutor
            from concurrent.futures.process import BrokenProcessPool
            context = multiprocessing.get_context("fork")
        except (ImportError, ValueError):
            why = "this platform cannot fork a process"
        else:
            try:
                with ProcessPoolExecutor(workers, mp_context=context, initializer=_one_thread) as pool:
                    results = pool.map(_run, range(len(jobs)), chunksize=1)
                    for row in results:
                        if isinstance(row, _Raised):
                            raised = True
                            results.close()                    # the jobs that no worker has taken are not run
                            break
                        rows.append(row)
            except BrokenProcessPool:
                why = "a worker process ended before it gave its row"
        if len(rows) < len(jobs):
            if raised:                                         # the job that raised in a worker is run here: if it raises again, that is what the caller gets
                rows.append(fn(jobs[len(rows)]))
                why = "a job raised in a worker and not in the calling process"      # it did not: the fault was the worker's alone, and the caller is told
            here = len(jobs) - len(rows)
            LOST.append({"why": why, "rows_computed_in_the_calling_process": here + (1 if raised else 0), "of": len(jobs)})
            warnings.warn(f"twinfit.pool: {why}. {here + (1 if raised else 0)} of {len(jobs)} rows are computed in the calling process, one "
                          "after another: the rows are the same and the time is longer.", RuntimeWarning, stacklevel=2)
            rows += [fn(j) for j in jobs[len(rows):]]
        return rows
    finally:
        _JOBS = []
