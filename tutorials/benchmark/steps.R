# Step recorder for the R side of the benchmark — the mirror of steps.py.
#
# Writes one JSON line per timed step. `as.numeric(Sys.time())` is epoch
# seconds, the same clock Python's `time.time()` reads, so the parent can align
# these boundaries with the RSS samples it took from outside the process.
#
# Deliberately hand-rolled JSON rather than jsonlite::toJSON: the log has to
# survive being written from a process that is about to run out of memory, and
# one `cat` per line with an open connection is the least that can go wrong.

step_log_open <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  env <- new.env(parent = emptyenv())
  env$fh <- file(path, open = "w")
  env
}

.step_emit <- function(log, name, t0, t1, anchor) {
  anchor_json <- if (is.null(anchor)) "null" else {
    if (is.character(anchor)) sprintf('"%s"', anchor) else
      format(anchor, scientific = FALSE, trim = TRUE)
  }
  cat(sprintf('{"step": "%s", "t0": %.6f, "t1": %.6f, "seconds": %.6f, "anchor": %s}\n',
              name, t0, t1, t1 - t0, anchor_json),
      file = log$fh)
  flush(log$fh)
}

# Time `expr`. The expression's value becomes the step's anchor if it is a
# single number or string, and is returned to the caller either way.
step <- function(log, name, expr) {
  t0 <- as.numeric(Sys.time())
  value <- force(expr)
  t1 <- as.numeric(Sys.time())
  anchor <- if (is.atomic(value) && length(value) == 1L &&
                (is.numeric(value) || is.character(value))) value else NULL
  .step_emit(log, name, t0, t1, anchor)
  invisible(value)
}

# Record a step timed by hand (library() calls, mostly).
step_mark <- function(log, name, seconds, anchor = NULL) {
  t1 <- as.numeric(Sys.time())
  .step_emit(log, name, t1 - seconds, t1, anchor)
}

step_log_close <- function(log) close(log$fh)
