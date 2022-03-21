[![Continous integration](https://github.com/Redmar-van-den-Berg/haploblock-shuffler/actions/workflows/ci.yml/badge.svg)](https://github.com/Redmar-van-den-Berg/haploblock-shuffler/actions/workflows/ci.yml)

# Haploblock-shuffler
Create all possible combinations of phased and unphased blocks in a vcf

------------------------------------------------------------------------
## Background
This tool takes a phase, unphased or partially phased VCF file, and generates
all possible combinations of phase blocks that are consistent with the phasing
that is present in the VCF file.


## Details
First, this tool reads all variants from a VCF file, and groups variants
together if they are compatible.
1. If a variant is phased (using the `PS` tag), it is only compatible with
   other phased variants that have the same phase ID.
2. Homozygous variants are always compatible with other variants, since they
   are part of every phase group
3. Heterozygous variants are only compatible when they are phased, and the
   phase ID matches.

To produce all possible combinations of grouped variants, haplotype-suffler
uses a counter to produce a binary pattern that determines which calls should
be modified. To modify a variant, we simply invert the order of the `GT` field,
so that `0/1` becomes `1/0`, or vice versa.

Since there are two alleles for every variant, we only have to produce half of
the possible VCF file, since the other half are mirror images (e.g. `0101` and
`1010`).

## Usage
```bash
haploblock-shuffler test.vcf output
```

To generate consensus fasta files from the output vcf files, bgzip and index
the output vcf files
```bash
cd output
for i in out_*.vcf; do
    bgzip $i
    tabix ${i}.gz
done
```
Then, generate the consensus using
```bash
samtools faidx $REFERENCE $REGION | bcftools consensus -H 1 out_0.vcf.gz > out_0_1.fa
samtools faidx $REFERENCE $REGION | bcftools consensus -H 2 out_0.vcf.gz > out_0_2.fa
```

## Limitations
This tool will generate `2^(n-1)` VCF files in the specified `output` folder,
where `n` is the number of phase blocks in the input VCF (see above). By
default, this is limited to 11 blocks, which means that at most 1024 files will
be created. This limit can be increased by using `--max-blocks`, but use with
caution.
