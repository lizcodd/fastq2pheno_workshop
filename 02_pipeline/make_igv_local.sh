#!/bin/bash
# Package the pipeline's coverage tracks and peaks for IGV on your own computer.
#
# The pipeline's own session (results/igv/narrow_peak/igv_session.xml) points at the server's copy of the genome,
# so it only opens on the server. This makes igv_local/ with the same kind of tracks and a session that uses
# IGV's built-in hg38 instead, and packs it into igv_local.tar (about 1.1 GB) to download.
#
# Run from 02_pipeline/ after the pipeline has finished (under a minute):  bash make_igv_local.sh
set -euo pipefail

out=igv_local
bwa=results/bwa
mkdir -p "$out"

# Links to the pipeline's files, under shorter names ("mRp" = both replicates merged, "mLb" = one replicate)
link() { ln -sfn "$(realpath "$1")" "$out/$2"; }
link $bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.bed  consensus_peaks.bed
for line in COLO320DM COLO320HSR; do
    link $bwa/merged_replicate/macs2/narrow_peak/$line.mRp.clN_peaks.narrowPeak    ${line}_peaks.narrowPeak
    link $bwa/merged_replicate/bigwig/$line.mRp.clN.bigWig                         $line.bigWig
    for rep in REP1 REP2; do
        link $bwa/merged_library/bigwig/${line}_$rep.mLb.clN.bigWig                ${line}_$rep.bigWig
    done
done

# The IGV session: genome, starting region, and one line per track (DM red, HSR blue).
# autoscaleGroup="1" puts all coverage tracks on the same y scale.
feature() { echo "    <Track clazz=\"org.broad.igv.track.FeatureTrack\" id=\"$1\" name=\"$2\" color=\"$3\" displayMode=\"COLLAPSED\" height=\"20\"/>"; }
coverage() { echo "    <Track clazz=\"org.broad.igv.track.DataSourceTrack\" id=\"$1\" name=\"$2\" color=\"$3\" renderer=\"BAR_CHART\" autoScale=\"true\" autoscaleGroup=\"1\" height=\"60\"/>"; }
{
    echo '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
    echo '<Session genome="hg38" locus="chr8:126300000-128200000" version="8">'
    echo '  <Resources>'
    for f in consensus_peaks.bed COLO320DM_peaks.narrowPeak COLO320HSR_peaks.narrowPeak COLO320DM.bigWig COLO320HSR.bigWig \
             COLO320DM_REP1.bigWig COLO320DM_REP2.bigWig COLO320HSR_REP1.bigWig COLO320HSR_REP2.bigWig; do
        echo "    <Resource path=\"$f\"/>"
    done
    echo '  </Resources>'
    echo '  <Panel name="DataPanel">'
    feature  consensus_peaks.bed          "consensus peaks"       "0,0,0"
    feature  COLO320DM_peaks.narrowPeak   "DM peaks"              "178,34,34"
    feature  COLO320HSR_peaks.narrowPeak  "HSR peaks"             "70,130,180"
    coverage COLO320DM.bigWig             "DM (both replicates)"  "178,34,34"
    coverage COLO320HSR.bigWig            "HSR (both replicates)" "70,130,180"
    coverage COLO320DM_REP1.bigWig        "DM rep 1"              "205,92,92"
    coverage COLO320DM_REP2.bigWig        "DM rep 2"              "205,92,92"
    coverage COLO320HSR_REP1.bigWig       "HSR rep 1"             "100,149,237"
    coverage COLO320HSR_REP2.bigWig       "HSR rep 2"             "100,149,237"
    echo '  </Panel>'
    echo '</Session>'
} > "$out/igv_session.xml"

cat > "$out/README.txt" <<'EOF'
COLO320DM vs COLO320HSR ATAC-seq tracks for IGV desktop (hg38).

1. Unpack this folder on your computer (keep the files together).
2. In IGV: File > Open Session... > igv_session.xml. It uses IGV's built-in hg38 genome, so no genome download.
3. It opens on the MYC amplicon (chr8). Type a gene or region in the search box to move.

Tracks: consensus peaks (all samples), each line's peaks, coverage of each line (both replicates merged) and of
each replicate. Coverage is nf-core/atacseq's bigWig, scaled to 1 million mapped reads. All coverage tracks share
one y scale. Heights are not corrected for copy number or signal-to-noise, so the amplicon towers over everything;
zoom to a normal region, or right-click > Set Data Range, to compare there.
EOF

# One file to download; -h stores the files the links point to, not the links
tar -chf igv_local.tar "$out"
ls -lh igv_local.tar
