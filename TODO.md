# TODO

## Docs

- Some non-updated docs, e.g. on the CLI:
   - "--verbose [silent|info|debug]” vs. "-v, --verbosity [info|debug]”
   - “hist  Plot (stacked) histograms.” vs. “viz  Plot histograms/densities (interactively)."

## General

- print the default values for CLI options when using -h (e.g. for `wgd ks -h` note the default value “fasttree” for option “-w, —wm”)

- report versions of third party software in logs

- add log file option

- log the total runtime

- use n_threads instead of n_cores

## `wgd blast`

- `wgd blast` sometimes doesn’t work correctly when “-n_threads" > 7 or so.
Seems like blastp crashes or doesn’t finish fully, but there are no errors
reported in the log output (just "INFO Blast done” comes immediately or very
quickly), and then mcl crashes or has very few entries. May be related to how
many cores and/or not enough memory requested.

- Check gene names and report issues

## `wgd ks`

- running `wgd ks` hangs or crashes when the path of the working directory is
  too long. More specifically codeml crashes (or silently hangs in the
background waiting for input) when the paths in “seqfile” or “outfile” are too
long. The symptoms for this can be different, e.g., it just hangs at "DEBUG
Codeml iteration 1 for […]” (because codeml is waiting for input), or it
crashes with “WARNING Codeml output file […].codeml not found”, etc.. Maybe
just replace the full path in “seqfile” and “outfile” with just the filename,
since it seems like codeml is being executed from within the temporary
directory already anyways. And maybe apply the same to fasttree and phyml.

- running `wgd ks` with phyml crashes. Looks like it’s crashing when it comes
  to a gene family with only two genes and tries to (unnecessarily) build the
tree. Skip running phyml for these.

- when running `wgd ks` occasionally got three “mv: Argument list too long” log outputs at the end after "INFO Moving files (preserve=True)”, and thus the codeml, msa, and tree folders are empty

- remove orphan rst rst1 rub files after wgd ks
