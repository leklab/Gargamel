version 1.0

workflow GargamelPipeline {

  input {
    File inputSamplesFile

    String clean_reads_py

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
        sample_name = sample[0],
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
        clean_reads_py = clean_reads_py
    }


  }

  output {
    Array[File] raw_bam_files = AlignReads.output_bam
    Array[File] bam_files = CleanReads.output_bam
    Array[File] bam_index = CleanReads.output_bam_index

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
    bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" ~{ref_fasta} \
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
  }

  command <<<
    python ~{clean_reads_py} -l 26 -i ~{raw_bam} -o - \
    | samtools view -Sh1 - \
    | samtools sort - > ~{sample_name}.bam

    samtools index ~{sample_name}.bam
  >>>

  runtime {
    cpus: 2
    requested_memory: 8000
  }

  output {
    File output_bam = "~{sample_name}.bam"
    File output_bam_index = "~{sample_name}.bam.bai"
  }
}

