"""End-to-end ``buildidx`` against a synthetic genome.

Small enough to run in seconds, but it exercises the whole path:
argument parsing, data-directory resolution, archive extraction, the
N/GC/repeat bit arrays, the null-seq index, and the FASTA index.
"""

import pytest


WINDOW = 100  # keep the (t+1)x(t+1) null-seq table small


def _build(run_cli, tiny_genome, cwd, extra=()):
    argv = ["buildidx", "-i", tiny_genome["archive"], "-g", tiny_genome["genome"],
            "-w", WINDOW, "-@", 1] + list(extra)
    return run_cli(argv, cwd=cwd)


def test_buildidx_from_within_the_data_dir(run_cli, tiny_genome, data_dir):
    """The README flow: ``cd data && gkmqc buildidx ...``.

    Under a non-editable install the old resolver sent this into
    site-packages and died in os.mkdir.
    """
    proc = _build(run_cli, tiny_genome, cwd=data_dir)
    assert proc.returncode == 0, proc.stdout

    genome = data_dir / tiny_genome["genome"]
    assert genome.is_dir(), "index was not built under the data directory"
    for chrom in tiny_genome["chroms"]:
        assert (genome / "fa" / ("%s.fa" % chrom)).is_file()
        for kind in ("cg", "na", "rp"):
            assert (genome / "bit" / ("%s.%s.bit" % (chrom, kind))).is_file()
        assert (genome / ("nidx_t%d" % WINDOW) / ("%s_pos.npy" % chrom)).is_file()
        assert (genome / ("nidx_t%d" % WINDOW) / ("%s_ptr.npz" % chrom)).is_file()


def test_buildidx_writes_fasta_index(run_cli, tiny_genome, data_dir):
    """The .fai is built up front, not left to the first evaluate run.

    Building it here means the directory is writable by construction and
    there is a single writer, so concurrent evaluate runs against a
    shared data directory cannot race to create it.
    """
    proc = _build(run_cli, tiny_genome, cwd=data_dir)
    assert proc.returncode == 0, proc.stdout

    fa_dir = data_dir / tiny_genome["genome"] / "fa"
    for chrom in tiny_genome["chroms"]:
        fai = fa_dir / ("%s.fa.fai" % chrom)
        assert fai.is_file(), "no .fai for %s" % chrom
        assert fai.read_text().startswith(chrom)


def test_buildidx_honours_explicit_data_dir(run_cli, tiny_genome, tmp_path):
    """-D wins, and buildidx may create the directory it is pointed at."""
    target = tmp_path / "not_created_yet"
    work = tmp_path / "somewhere_else"
    work.mkdir()
    proc = _build(run_cli, tiny_genome, cwd=work, extra=["-D", target])
    assert proc.returncode == 0, proc.stdout
    assert (target / tiny_genome["genome"] / "fa").is_dir()


def test_buildidx_honours_env_var(run_cli, tiny_genome, tmp_path):
    target = tmp_path / "via_env"
    target.mkdir()
    work = tmp_path / "elsewhere"
    work.mkdir()
    proc = _build(run_cli, tiny_genome, cwd=work,
                  extra=[])  # no -D; rely on the environment
    # re-run with the env var set
    proc = run_cli(
        ["buildidx", "-i", tiny_genome["archive"], "-g", tiny_genome["genome"],
         "-w", WINDOW, "-@", 1],
        cwd=work, env={"GKMQC_DATA_DIR": str(target)},
    )
    assert proc.returncode == 0, proc.stdout
    assert (target / tiny_genome["genome"] / "fa").is_dir()


def test_buildidx_is_rerunnable(run_cli, tiny_genome, data_dir):
    """A second run must not fail on directories that already exist."""
    first = _build(run_cli, tiny_genome, cwd=data_dir)
    assert first.returncode == 0, first.stdout
    second = _build(run_cli, tiny_genome, cwd=data_dir)
    assert second.returncode == 0, second.stdout
