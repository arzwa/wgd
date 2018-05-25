# TODO

## Docs

## General

## `wgd blast`

- `wgd blast` sometimes doesn’t work correctly when “-n_threads" > 7 or
so. Seems like blastp crashes or doesn’t finish fully, but there are no
errors reported in the log output (just "INFO Blast done” comes
immediately or very quickly), and then mcl crashes or has very few
entries. May be related to how many cores and/or not enough memory
requested.

- Check gene names and report issues

## `wgd ks`

- running `wgd ks` with phyml crashes. Looks like it’s crashing when it
comes to a gene family with only two genes and tries to (unnecessarily)
build the tree. Skip running phyml for these.

- Add option to use fixed kappa (as in Lynch & Conery 2003)

