# Per-variant VAF vs genotype concordance for mutsig VAF-selected VCFs.
# Input: bcftools query lines with columns:
#   CHROM POS REF ALT GT F R POSITION
BEGIN {
    FS = "\t"
    OFS = "\t"
    print "VariantIt", "VAF", "GT", "Verdict"
    d = 0
    n = 0
}
{
    vid = $1 ":" $2 ":" $3 ":" $4
    f = $6 + 0
    r = $7 + 0
    p = $8 + 0
    if (p <= 0 || $8 == ".")
        vaf = "NA"
    else
        vaf = (f + r) / p
    gt = $5
    if (gt == "0/0" || gt == "0|0") {
        ver = "Discordant"
        d++
    } else
        ver = "Concordant"
    n++
    print vid, vaf, gt, ver
}
END {
    if (n > 0)
        print "# discordant_proportion", d / n, d, n
    else
        print "# discordant_proportion", "NA", "NA", "NA"
}
