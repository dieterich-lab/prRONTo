rule samtools_preprocess_bam:
  input: bam=join_path("data/bams/raw/{filename}.bam"),
         bai=join_path("data/bams/raw/{filename}.bam.bai"),
         ref=REF_FASTA,
         regions=REGIONS,
  output: join_path("results/data/bams/preprocessed/{filename}.bam"),
  params:
    filter=config["preprocess"]["filter"],
    calmd=config["preprocess"]["calmd"],
  log: join_path("logs/samtools/preprocess/{filename}.log")
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
  input: join_path("results/data/bams/preprocessed/{bam}")
  output: join_path("results/downsampling/bams/seed~{seed}_reads~{reads}/{bam}")
  log: join_path("logs/samtools/downsample_bam/seed~{seed}_reads~{reads}/{bam}.log")
  script: "scripts/sample_bam.py"


rule samtools_mix_bams:
  input: condA=join_path("results/downsampling/bams/seed~{seed}_reads~{reads}/{bamA}"),
         condB=join_path("results/downsampling/bams/seed~{seed}_reads~{reads}/{bamB}"),
  output: join_path("results/mixing/bams/seed~{seed}_reads~{reads}/cond1~{fractionA}_and_cond2~{fractionB}.bam")
  log: join_path("logs/samtools/mix_bams/seed~{seed}_reads~{reads}/cond1~{fractionA}_and_cond2~{fractionB}.log")
  params:
    condA=lambda wildcards: f"--subsample-seed {wildcards.seed} --subsample {wildcards.fractionA}",
    condB=lambda wildcards: f"--subsample-seed {wildcards.seed} --subsample {wildcards.fractionB}",
  shell: """
    samtools merge -c -p \
        {output} \
        <(samtools view {params.condA} {input.condA}) \
        <(samtools view {params.condB} {input.condB}) \
        2> {log}
    """


rule samtools_index_bam:
  input: "{filename}.bam"
  output: "{filename}.bam.bai"
  shell: "samtools index {input}"
