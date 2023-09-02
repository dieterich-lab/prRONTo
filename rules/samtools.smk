rule samtools_preprocess_bam:
  input: bam=join_path("data/bams/{bam}"),
         ref=REF_FASTA,
         regions=REGIONS,
  output: join_path("results/data/preprocessed/bams/{bam}"),
  params:
    filter=config["preprocess"]["filter"],
    calmd=config["preprocess"]["calmd"],
  log: join_path("logs/samtools/preprocess/{bam}")
  run:
    filter_cmd = "samtools view {params.filter} -b {input}"
    calmd_cmd = "samtools calmd -b {input.ref} > {output}"

    cmds = [filter_cmd, ]
    if params.callmd:
      cmds.append(calmd_cmd)
    cmd = "|".join(cmds)

    shell("( " + cmd " > {output} ) 2> {log}")


rule samtools_downsample_bam:
  input: join_path("results/data/preprocessed-bams/{bam}")
  output: join_path("results/downsampling/bams/seed~{seed}_reads~{reads}/{bam}")
  log: join_path("logs/samtools/downsample_bam/seed~{seed}_reads~{reads}/{bam}")
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
  input: "{bam}"
  output: "{bam}.bai"
  shell: """
    samtools index {input}
  """
