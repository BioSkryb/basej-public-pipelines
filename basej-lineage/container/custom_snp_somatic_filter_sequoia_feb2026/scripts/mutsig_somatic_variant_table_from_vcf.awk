# Variant table (example.txt-style) from bcftools query: CHROM POS REF ALT
# Classifies each allele as SNP, DBS, or INDEL and emits one row per allele.
#   SNP:   len(REF)==1 && len(ALT)==1  → mut_type="SNP",   pos_end=pos_start
#   DBS:   len(REF)==2 && len(ALT)==2  → mut_type="DBS",   pos_end=pos_start
#   INDEL: otherwise                   → mut_type="INDEL", pos_end=pos_start+len(REF)-1
# Multi-allelic ALT fields (comma-separated) are split and classified independently.
# Usage: bcftools query ... | awk -v project=... -v sample=... -v genome=... -f thisfile.awk
BEGIN {
    FS = "\t"
    OFS = "\t"
    print "Project", "Sample", "ID", "Genome", "mut_type", "chrom", "pos_start", "pos_end", "ref", "alt", "Type"
}
{
    chrom = $1
    pos = $2
    ref = $3
    alts = $4
    n = split(alts, a, ",")
    for (i = 1; i <= n; i++) {
        alt = a[i]
        if (length(ref) == 1 && length(alt) == 1) {
            mut_type = "SNP"
            pos_end = pos
        } else if (length(ref) == 2 && length(alt) == 2) {
            mut_type = "DBS"
            pos_end = pos
        } else {
            mut_type = "INDEL"
            pos_end = pos + length(ref) - 1
        }
        print project, sample, ".", genome, mut_type, chrom, pos, pos_end, ref, alt, "SOMATIC"
    }
}
