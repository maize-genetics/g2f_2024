#!/usr/bin/env bash
# Confirm that the original sequencing files listed in the raw-data manifest
# still exist, at the manifest size, on the source storage host. Read-only:
# one SSH session running `stat` (you may be asked for your password).
#
#   bash scripts/check_remote_copies.sh MANIFEST USER@HOST [OUT]
#
# Writes OUT (default ../data/remote_copies.tsv; set it as remote_copy_table
# in the top-level configuration): source_path, exists, size, mtime, checked_at.
# Then refresh the release report:
#   python3 scripts/g2f_tool.py raw-release --config <top-level config>
set -euo pipefail
cd "$(dirname "$0")/.."
[[ $# -ge 2 ]] || { echo "usage: $0 MANIFEST USER@HOST [OUT]" >&2; exit 2; }
MANIFEST=$1
HOST=$2
OUT=${3:-../data/remote_copies.tsv}
[[ -f "$MANIFEST" ]] || { echo "manifest not found: $MANIFEST" >&2; exit 1; }
paths=$(awk -F'\t' 'NR>1 && $6 ~ /^\// {print $6}' "$MANIFEST" | tr -d '\r')
[[ -n "$paths" ]] || { echo "no source paths in $MANIFEST" >&2; exit 1; }
now=$(date '+%Y-%m-%dT%H:%M:%S')
tmp=$(mktemp "${OUT}.XXXX")
printf 'source_path\texists\tsize\tmtime\tchecked_at\n' > "$tmp"
# One remote shell; paths are passed on stdin (no quoting issues).
printf '%s\n' $paths | ssh -o ConnectTimeout=20 "$HOST" '
    while IFS= read -r p; do
        if s=$(stat -c "%s %Y" -- "$p" 2>/dev/null); then
            printf "%s\tyes\t%s\t%s\n" "$p" "${s% *}" "${s#* }"
        else
            printf "%s\tno\t\t\n" "$p"
        fi
    done' | awk -F'\t' -v t="$now" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,t}' >> "$tmp"
n=$(($(wc -l < "$tmp") - 1)); want=$(printf '%s\n' $paths | wc -l)
[[ $n -eq $want ]] || { echo "expected $want rows, got $n; not written" >&2; rm -f "$tmp"; exit 1; }
mv "$tmp" "$OUT"
echo "wrote $OUT: $(awk -F'\t' 'NR>1 && $2=="yes"' "$OUT" | wc -l) of $n present on $HOST"
# Size check against the manifest
awk -F'\t' 'NR==FNR {if (FNR>1) m[$6]=$7; next} FNR>1 && ($2!="yes" || $3!=m[$1]) {print "  MISMATCH/MISSING:", $1, $2, $3, "(manifest", m[$1] ")"}' \
    <(tr -d '\r' < "$MANIFEST") "$OUT"
