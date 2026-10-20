#!/bin/bash
# Thread-scaling benchmark: cpenc vs zstd/xz/gzip on the remote ARM box.
# Emits pipe-delimited rows: FILE|TOOL|THREADS|RATIO|CSIZE|CMBps|DMBps
cd "$HOME/bench" || exit 1
TMP="$HOME/btmp"; mkdir -p "$TMP"
OUT="$TMP/c.bin"

NOW=0
now() { local t=${EPOCHREALTIME}; NOW=$(( ${t%.*}*1000000 + 10#${t#*.} )); }

# best-of-2 compress: $1 = command string, stdout captured to $OUT. echoes best us
runbest_c() {
  local best=999999999999 i s el
  for i in 1 2; do
    now; s=$NOW
    eval "$1" > "$OUT"
    now; el=$(( NOW - s ))
    (( el < best )) && best=$el
  done
  echo "$best"
}

# best-of-2 decompress: $1 = command string reading $OUT, output discarded
runbest_d() {
  local best=999999999999 i s el
  for i in 1 2; do
    now; s=$NOW
    eval "$1" > /dev/null
    now; el=$(( NOW - s ))
    (( el < best )) && best=$el
  done
  echo "$best"
}

emit() { # file n tool threads csize cus dus
  awk -v f="$1" -v n="$2" -v t="$3" -v th="$4" -v sz="$5" -v c="$6" -v d="$7" \
    'BEGIN{printf "%s|%s|%s|%.4f|%d|%.1f|%.1f\n", f,t,th,n/sz,sz,n/c,n/d}'
}

# $1=file $2=n $3=tool $4=threads $5=compress-cmd $6=decompress-cmd
stream() {
  local cus sz dus
  cus=$(runbest_c "$5")
  sz=$(stat -c %s "$OUT")
  dus=$(runbest_d "$6")
  emit "$1" "$2" "$3" "$4" "$sz" "$cus" "$dus"
}

echo "### start $(date -u +%FT%TZ)"

for f in big_rec.bin big_text.bin tab8.bin big_rand.bin; do
  n=$(stat -c %s "$f")
  echo "### file $f n=$n"

  # --- cpenc (best of 2 invocations, each internally best-of-3) ---
  for th in 1 2 4; do
    r1=$(CPENC_THREADS=$th ./spd_new "$f" 3)
    r2=$(CPENC_THREADS=$th ./spd_new "$f" 3)
    echo "$r1|$r2" | awk -F'|' -v f="$f" -v th="$th" '{
      split($1,a," "); split($2,b," ");
      gsub("x","",a[1]); gsub("x","",b[1]);
      r=(a[1]>b[1])?a[1]:b[1]; c=(a[2]>b[2])?a[2]:b[2]; d=(a[3]>b[3])?a[3]:b[3];
      printf "%s|cpenc|%s|%.4f|%d|%.1f|%.1f\n", f, th, r, 0, c, d }'
  done

  # --- zstd -19 ---
  for th in 1 2 4; do
    stream "$f" "$n" "zstd-19" "$th" "zstd -19 -T$th -c $f" "zstd -d -c $OUT"
  done

  # --- zstd -3 ---
  for th in 1 2 4; do
    stream "$f" "$n" "zstd-3" "$th" "zstd -3 -T$th -c $f" "zstd -d -c $OUT"
  done

  # --- xz -6 (default single block for these sizes) ---
  for th in 1 2 4; do
    stream "$f" "$n" "xz-6" "$th" "xz -6 -T$th -c $f" "xz -d -c $OUT"
  done
  # xz -T0 (all 4 cores)
  stream "$f" "$n" "xz-6" "0(all)" "xz -6 -T0 -c $f" "xz -d -c $OUT"
  # xz -6 -T4 with explicit 1MiB blocks (forces multi-block -> -T can help)
  stream "$f" "$n" "xz-6-1Mblk" "4" "xz -6 -T4 --block-size=1MiB -c $f" "xz -d -c $OUT"
  # block count of the -T1 vs -T4 vs 1MiB-block outputs
  xz -6 -T1 -c "$f" > "$TMP/b1.xz"
  xz -6 -T4 -c "$f" > "$TMP/b4.xz"
  xz -6 -T4 --block-size=1MiB -c "$f" > "$TMP/b4m.xz"
  bl1=$(xz -l --robot "$TMP/b1.xz" | awk -F'\t' '$1=="file"{print $3}')
  bl4=$(xz -l --robot "$TMP/b4.xz" | awk -F'\t' '$1=="file"{print $3}')
  bl4m=$(xz -l --robot "$TMP/b4m.xz" | awk -F'\t' '$1=="file"{print $3}')
  echo "### xzblocks $f T1=$bl1 T4=$bl4 T4-1MiB=$bl4m"
  rm -f "$TMP"/b1.xz "$TMP"/b4.xz "$TMP"/b4m.xz

  # --- gzip -6 (single-threaded baseline) ---
  stream "$f" "$n" "gzip-6" "1" "gzip -6 -c $f" "gzip -d -c $OUT"

  # --- integrity check: decompress each produced stream and cmp ---
  for spec in "zstd -19 -T1|zstd -d" "zstd -3 -T1|zstd -d" "xz -6 -T1|xz -d" "gzip -6|gzip -d"; do
    e=${spec%%|*}; d=${spec##*|}
    $e -c "$f" > "$OUT"
    if $d -c "$OUT" | cmp -s - "$f"; then v=OK; else v=MISMATCH; fi
    echo "### verify $f ${e%% *} $v"
  done
done

echo "### done $(date -u +%FT%TZ)"
rm -f "$OUT" "$TMP"/b*.xz
rmdir "$TMP" 2>/dev/null
