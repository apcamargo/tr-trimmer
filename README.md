# tr-trimmer

Identify and trim terminal repeats from sequences in FASTA files.

## Usage

```
tr-trimmer [OPTIONS] [INPUT]...
```

### Arguments

| Argument | Description |
|:---------|:------------|
| `[INPUT]...` | Input file(s). Use `-` for stdin [default: `-`] |

### Options

| Option | Description |
|:-------|:------------|
| `-h`, `--help` | Print help |
| `-V`, `--version` | Print version |

### Terminal repeat identification

| Option | Description |
|:-------|:------------|
| `-i`, `--enable-itr-identification` | Identify inverted terminal repeats (ITRs) from sequences. |
| `-d`, `--disable-dtr-trimming` | Disable identification of direct terminal repeats (DTRs). Requires `--enable-itr-identification`. |
| `-l`, `--min-length <MIN_LENGTH>` | Minimum length of terminal repeat. Default: `21`. |

### Terminal repeat filtering

| Option | Description |
|:-------|:------------|
| `-c`, `--ignore-low-complexity` | Ignore terminal repeats with a high proportion of low-complexity sequences. |
| `--max-low-complexity-frac <MAX_LOW_COMPLEXITY_FRAC>` | Maximum fraction of the terminal repeat length comprised of low-complexity sequences. Default: `0.5`. |
| `-n`, `--ignore-ambiguous` | Ignore terminal repeats with a high proportion of ambiguous bases (e.g., `N`). |
| `--max-ambiguous-frac <MAX_AMBIGUOUS_FRAC>` | Maximum fraction of the terminal repeat length comprised of ambiguous bases. Default: `0.5`. |

### Output

| Option | Description |
|:-------|:------------|
| `-x`, `--exclude-non-tr-seqs` | Retain only sequences for which terminal repeats were identified. |
| `-a`, `--include-tr-info` | Add terminal repeat information to sequence headers (e.g., `tr=dtr tr_length=55`). |
| `-t`, `--disable-trimming` | Disable trimming of terminal repeats. Can be used with `--include-tr-info` or `--exclude-non-tr-seqs` to identify and report sequences with terminal repeats without modifying the sequences. |
