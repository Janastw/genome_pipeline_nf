process EXTRACT_16S_RRNA {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fna), path(gff)

    output:
    tuple val(meta), path("${meta.id}-16S.fna"), emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    gfffile=${gff}
    fastafile=${fna}

    grep '16S' \$gfffile > 16S-gff.gff
    bedtools getfasta -fi \$fastafile -bed 16S-gff.gff -fo 16S-fasta.fna
    grep -m 1 ">" 16S-fasta.fna | sed 's/>//g' > 16S-id.txt
    xargs samtools faidx 16S-fasta.fna < 16S-id.txt > ${meta.id}-16S.fna

    rm 16S-gff.gff
    rm 16S-fasta.fna
    rm 16S-fasta.fna.fai
    rm 16S-id.txt

    cat > versions.yml <<- 'END_VERSIONS'
        "${task.process}":
            bedtools: \$(bedtools --version | sed 's/bedtools v//g')
            samtools: \$(samtools --version | head -1 | sed 's/samtools //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}-16S.fna versions.yml
    """
}
