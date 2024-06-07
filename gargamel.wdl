version 1.0

workflow GargamelPipeline {

  input {
    String gene
    File inputSamplesFile

    String clean_reads_py
    String site_metrics_py
    String count_wt_blocks_py

    File ref_fasta
    File ref_fasta_index

    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

  }

  Array[Array[String]] inputSamples = read_tsv(inputSamplesFile)

  scatter (sample in inputSamples) {
    call AlignReads {
      input:
        sample_name = gene + '_' + sample[0],
        fastq1 = sample[1],
        fastq2 = sample[2],
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }

    call CleanReads {
      input:
        sample_name = sample[0],
        raw_bam = AlignReads.output_bam,
        clean_reads_py = clean_reads_py,
        overlap_len = sample[5]
    }

    call VariantCounts {
      input:
        sample_name = gene + '_' + sample[0],
        gene = gene,
        input_bam = CleanReads.output_bam,
        input_bam_index = CleanReads.output_bam_index,
        site_metrics_py = site_metrics_py,
        start_block = sample[3],
        end_block = sample[4]
    }

    call CountWTBlocks {
      input:
        sample_name = sample[0],
        input_bam = CleanReads.output_bam,
        input_bam_index = CleanReads.output_bam_index,
        count_wt_blocks_py = count_wt_blocks_py
    }

  }

  output {
    Array[File] raw_bam_files = AlignReads.output_bam
    Array[File] bam_files = CleanReads.output_bam
    Array[File] bam_index = CleanReads.output_bam_index
    Array[File] block_variant_counts = VariantCounts.block_variant_counts
    Array[File] wt_block_counts = CountWTBlocks.wt_block_counts

  }

}

task AlignReads {
  input {
    String fastq1
    String fastq2
    String sample_name

    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    File ref_fasta
    File ref_fasta_index

  }

  command <<<
    bwa mem -M -t 8 -R "@RG\tID:~{sample_name}\tPL:ILLUMINA\tSM:~{sample_name}" ~{ref_fasta} \
    ~{fastq1} ~{fastq2} 2> >(tee ~{sample_name}.bwa.stderr.log >&2) \
    | samtools view -Sh1 -F 0x900 \
    | samtools sort -n - > ~{sample_name}_raw.bam
  >>>

  runtime {
    cpus: 8
    requested_memory: 16000
  }

  output {
    File output_bam = "~{sample_name}_raw.bam"
  }
}


task CleanReads {
  input {  
    File raw_bam
    String clean_reads_py
    String sample_name
    String overlap_len
  }

  command <<<
    python ~{clean_reads_py} -l ~{overlap_len} -i ~{raw_bam} -o - \
    | samtools view -Sh1 - \
    | samtools sort - > ~{sample_name}.bam

    samtools index ~{sample_name}.bam
  >>>

  runtime {
    cpus: 1
    requested_memory: 8000
  }

  output {
    File output_bam = "~{sample_name}.bam"
    File output_bam_index = "~{sample_name}.bam.bai"
  }
}

task VariantCounts {
  input {  
    String sample_name
    String gene
    File input_bam
    File input_bam_index
    String site_metrics_py
    String start_block
    String end_block
  }

  command <<<
    python ~{site_metrics_py} -i ~{input_bam} -c ~{gene} -s ~{start_block} -e ~{end_block} -q 37 > ~{sample_name}_variant_counts.tsv
  >>>

  runtime {
    cpus: 1
    requested_memory: 64000
  }

  output {
    File block_variant_counts = "~{sample_name}_variant_counts.tsv"
  }
}

task CountWTBlocks {
  input {  
    String sample_name
    File input_bam
    File input_bam_index
    String count_wt_blocks_py
  }

  command <<<
    samtools sort -n -O SAM ~{input_bam} | python ~{count_wt_blocks_py} -i - -s ~{sample_name} > ~{sample_name}_wt_blocks.tsv
  >>>

  runtime {
    cpus: 1
    requested_memory: 16000
  }

  output {
    File wt_block_counts = "~{sample_name}_wt_blocks.tsv"
  }
}

