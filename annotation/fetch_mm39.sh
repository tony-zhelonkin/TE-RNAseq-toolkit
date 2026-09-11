#!/usr/bin/env bash
# fetch_mm39.sh — Layer 0 reference sources for mouse GRCm39 / GENCODE vM37.
#
# Writes an immutable dated snapshot and flips `current` at it.  Runs at either site:
#   workstation   --dest /data2/users/shared/refcache/mm39
#   CRI Randi     --dest /gpfs/data/rathmell-lab/data/refdata/mouse/GRCm39
#
# Every file is md5-verified against a pin recorded below.  A mismatch is a hard error:
# the whole TE pipeline is coordinate-exact, so a silently changed source is worse than
# a failed fetch.
set -euo pipefail

DEST=""
DATE=""
TE_GTF_SRC="${TE_GTF_SRC:-}"
DRY=0

usage() {
  cat <<'USAGE'
usage: fetch_mm39.sh --dest DIR [--date YYYYMMDD] [--te-gtf PATH] [--dry-run]

  --dest     parent directory; the snapshot lands at DIR/mm39_<date>/
  --date     snapshot tag date (default: today, UTC)
  --te-gtf   path to GRCm39_Ensembl_rmsk_TE.gtf.gz (see "TE GTF" in README.md)
  --dry-run  print what would happen

The GENCODE files download from a stable public URL.  The TE GTF does not — see README.md.
USAGE
}

while [ $# -gt 0 ]; do
  case "$1" in
    --dest)    DEST="$2"; shift 2 ;;
    --date)    DATE="$2"; shift 2 ;;
    --te-gtf)  TE_GTF_SRC="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
[ -n "$DEST" ] || { echo "--dest is required" >&2; usage; exit 2; }
DATE="${DATE:-$(date -u +%Y%m%d)}"

GENCODE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37

# name -> url
declare -A URL=(
  [GRCm39.genome.fa.gz]="$GENCODE/GRCm39.genome.fa.gz"
  [gencode.vM37.primary_assembly.annotation.gtf.gz]="$GENCODE/gencode.vM37.primary_assembly.annotation.gtf.gz"
)
# name -> md5.  Pinned from the copies that produced every prior run (13036-DM, 14839-DM).
declare -A MD5=(
  [GRCm39.genome.fa.gz]=c4dfe32c560ffabf5068f57449a14511
  [gencode.vM37.primary_assembly.annotation.gtf.gz]=702e27a8727db66cd23e9f66f5aa1d38
  [GRCm39_Ensembl_rmsk_TE.gtf.gz]=1d6770be131fc1b47402a633f4401cbd
)

SNAP="$DEST/mm39_$DATE"
echo "snapshot: $SNAP"
[ "$DRY" -eq 1 ] && { printf '  would fetch %s\n' "${!URL[@]}"; echo "  would verify and place GRCm39_Ensembl_rmsk_TE.gtf.gz from ${TE_GTF_SRC:-<--te-gtf required>}"; exit 0; }

mkdir -p "$SNAP"
cd "$SNAP"

verify() {  # verify <file>
  local f="$1" want="${MD5[$1]}" got
  got=$(md5sum "$f" | cut -d' ' -f1)
  [ "$got" = "$want" ] || { echo "MD5 MISMATCH $f: got $got want $want" >&2; exit 1; }
  echo "  ok  $f  $got"
}

for f in "${!URL[@]}"; do
  [ -f "$f" ] || curl -fsSL --retry 3 -o "$f" "${URL[$f]}"
  verify "$f"
done

TE=GRCm39_Ensembl_rmsk_TE.gtf.gz
if [ ! -f "$TE" ]; then
  [ -n "$TE_GTF_SRC" ] || { echo "ERROR: $TE absent and --te-gtf not given. See README.md." >&2; exit 1; }
  cp -n "$TE_GTF_SRC" "$TE"
fi
verify "$TE"

cat > MANIFEST.json <<JSON
{
  "source": "mm39",
  "snapshot_tag": "mm39_$DATE",
  "downloaded_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "assembly": "GRCm39",
  "gene_annotation": "GENCODE vM37",
  "te_annotation": "TEtranscripts mm39_rmsk (Ensembl contig names)",
  "contig_namespace": {
    "GRCm39.genome.fa.gz": "chr-prefixed",
    "gencode.vM37.primary_assembly.annotation.gtf.gz": "chr-prefixed",
    "GRCm39_Ensembl_rmsk_TE.gtf.gz": "unprefixed"
  },
  "upstream": {
    "gencode": "$GENCODE",
    "te_gtf": "TEtranscripts curated TE_GTF, distributed via Dropbox from https://www.mghlab.org/software/tetranscripts (no stable per-file URL)"
  },
  "md5": {
$(for f in GRCm39.genome.fa.gz gencode.vM37.primary_assembly.annotation.gtf.gz "$TE"; do
    printf '    "%s": "%s",\n' "$f" "${MD5[$f]}"; done | sed '$ s/,$//')
  },
  "file_count": 3,
  "total_size_bytes": $(du -bc GRCm39.genome.fa.gz gencode.vM37.primary_assembly.annotation.gtf.gz "$TE" | tail -1 | cut -f1)
}
JSON

ln -sfn "mm39_$DATE" "$DEST/current"
echo "current -> mm39_$DATE"
cat MANIFEST.json
