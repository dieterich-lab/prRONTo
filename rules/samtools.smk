rule samtools_preprocess_bam:
  input: bam=join_path("data/bams/raw/{filename}.bam"),
         bai=join_path("data/bams/raw/{filename}.bam.bai"),
         ref=REF_FASTA,
         regions=REGIONS,
  output: join_path("results/data/bams/preprocessed/{filename}.bam"),
  params:
    filter=config["preprocess"]["filter"],
    calmd=config["preprocess"]["calmd"],
  log: join_path("logs/samtools/preprocess_bam/{filename}.log") # TODO - where to put is in run
  run:
    with open(input.regions, "r") as f:
      regions = [line.strip() for line in f.readlines()]

    filter_cmd = "samtools view {params.filter} -b {input.bam} " + " ".join(regions)
    calmd_cmd = "samtools calmd -b /dev/stdin {input.ref}"

    cmds = [filter_cmd, ]
    if params.calmd:
      cmds.append(calmd_cmd)
    cmd = " | ".join(cmds)

    shell("( " + cmd + " > {output} ) 2> {log}")


rule samtools_downsample_bam:
  input: bam=join_path("results/data/bams/preprocessed/{filename}.bam"),
         coverage=join_path("results/data/bams/preprocessed/{filename}_coverage.tsv"),
  output: join_path("results/downsampling/bams/seed~{seed}_reads~{reads}/{filename}.bam")
  log: join_path("logs/samtools/downsample_bam/seed~{seed}_reads~{reads}/{filename}.log") # TODO - where to put is in run
  run:
    df = pd.read_csv(input.coverage, sep = "\t")
    read_count = sum(df["numreads"])
    fraction = round(int(wildcards.reads) / read_count, 3)
    if fraction > 1:
      # TODO up sample
      breakpoint()
    cmd = f"samtools view --subsample {fraction} --subsample-seed {wildcards.seed} " + "-o {output} {input.bam}"
    shell(cmd + " 2> {log}")


# TODO for mixing
#rule samtools_mix_bams:
#  input: condA=join_path("results/downsampling/bams/seed~{seed}_reads~{reads}/{bamA}"),
#         condB=join_path("results/downsampling/bams/seed~{seed}_reads~{reads}/{bamB}"),
#  output: join_path("results/mixing/bams/seed~{seed}_reads~{reads}/cond1~{fractionA}_and_cond2~{fractionB}.bam"),
#  log: join_path("logs/samtools/mix_bams/seed~{seed}_reads~{reads}/cond1~{fractionA}_and_cond2~{fractionB}.log"),
#  params:
#    condA=lambda wildcards: f"--subsample-seed {wildcards.seed} --subsample {wildcards.fractionA}",
#    condB=lambda wildcards: f"--subsample-seed {wildcards.seed} --subsample {wildcards.fractionB}",
#  shell: """
#    samtools merge -c -p \
#        {output} \
#        <(samtools view {params.condA} {input.condA}) \
#        <(samtools view {params.condB} {input.condB}) \
#        2> {log}
#    """


rule samtools_index_bam:
  input: "{filename}.bam"
  output: "{filename}.bam.bai"
  shell: "samtools index {input}"


rule samtools_coverage:
  input: "{filename}.bam",
  output: "{filename}_coverage.tsv"
  shell: "samtools coverage {input} > {output}"


rule samtools_stats:
  input: bam="{filename}.bam",
         ref=REF_FASTA,
  output: "{filename}_stats.tsv"
  shell: "samtools stats --ref-seq {input.ref} {input.bam} > {output}"


def raw_fnames(suffix):
  fnames = []
  for condition in SAMPLES["condition"].unique().tolist():
    for replicate in range(1, get_replicates(condition) + 1):
      fnames.append(join_path(f"data/bams/raw/cond{condition}_rep{replicate}{suffix}"))

  return fnames


def preprocessed_fnames(suffix):
  fnames = []
  for condition in SAMPLES["condition"].unique().tolist():
    for replicate in range(1, get_replicates(condition) + 1):
      fnames.append(join_path(f"results/data/bams/preprocessed/cond{condition}_rep{replicate}{suffix}"))

  return fnames


rule samtools_stats_extract_section:
  input: join_path("{prefix}/cond{condition}_rep{replicate}_stats.tsv"),
  output: join_path("{prefix}/cond{condition}_rep{replicate}_stats_{section}.tsv"),
  shell: "cat {input} | grep ^{wildcards.section} | cut -f 2- > {output}"


def merge_sections(input, cols, output):
  def helper(fname):
    df = pd.read_csv(fname, sep="\t", names = cols)
    df["fname"] = fname
    result = re.search(r".+/(raw|preprocessed)/cond([12])_rep([0-9]+)_stats_.+.tsv", fname)
    if not result:
      breakpoint()
    bam_type, cond, repl = result.groups()
    df["bam_type"] = bam_type
    df["condition"] = cond
    df["replicate"] = repl
    return df

  dfs = [helper(fname) for fname in input]
  df = pd.concat(dfs, ignore_index=True)
  df.to_csv(output[0], sep="\t", index=False)


rule merge_stats_section:
  input: raw_fnames("_stats_{section}.tsv") + preprocessed_fnames("_stats_{section}.tsv"),
  output: join_path("results/samtools/stats/merged_{section}.tsv"),
  params:
    columns=lambda wildcards: STATS_SECTION2COLUMNS[wildcards.section]
  run:
    merge_sections(input, params.columns, output)


rule merge_stats_sn:
  input: raw_fnames("_stats_SN.tsv") + preprocessed_fnames("_stats_SN.tsv"),
  output: join_path("results/samtools/stats/merged_SN.tsv"),
  run:
    def helper(fname):
      df = pd.read_csv(fname, sep="\t", names=["key", "value", "comment"])
      df["key"] = df["key"].str.replace(":", "")
      df = df[["key", "value"]].set_index("key", drop=True).transpose()
      df["fname"] = fname
      result = re.search(r".+/(raw|preprocessed)/cond([12])_rep([0-9]+)_stats_.+.tsv", fname)
      if not result:
        breakpoint()
      bam_type, cond, repl = result.groups()
      df["bam_type"] = bam_type
      df["condition"] = cond
      df["replicate"] = repl

      return df

    dfs = [helper(fname) for fname in input]
    df = pd.concat(dfs, ignore_index=True)
    df = df.convert_dtypes()
    df.to_csv(output[0], sep="\t", index=False)


def all_bams():
  fnames = raw_fnames("_coverage.tsv") + preprocessed_fnames("_coverage.tsv")

  for condition in SAMPLES["condition"].unique().tolist():
    for replicate in range(1, get_replicates(condition) + 1):
      if "downsampling" in config:
        for seed in config["downsampling"]["seed"]:
          for reads in config["downsampling"]["reads"]:
            fnames.append(
              join_path("results/downsampling/bams",
                        f"seed~{seed}_reads~{reads}",
                        f"cond{condition}_rep{replicate}_coverage.tsv"))

  return fnames


# TODO rename to coverage
rule samtools_read_summary:
  input: all_bams()
  output: join_path("results/read_summary.tsv"),
  log: join_path("logs/samtools/read_summary.log"), # TODO - where to put in run
  run:
    dfs = []
    for fname in input:
      df = pd.read_csv(fname, sep="\t")
      df["fname"] = fname
      result = re.search(r"bams/([^/]+)/cond(\d+)_rep(\d+)_coverage.tsv$", fname)
      parameters, condition, replicate = result.groups()
      df["parameters"] = parameters
      df["condition"] = condition
      df["replicate"] = replicate
      dfs.append(df)
    df = pd.concat(dfs)
    df.to_csv(output[0], sep="\t", index=False)
