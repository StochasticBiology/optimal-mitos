#!/usr/bin/env bash
set -euo pipefail

gif="$1"
avi="$2"
out="$3"

FPS=30

tmpdir=$(mktemp -d)
gif_dir="$tmpdir/gif"
avi_dir="$tmpdir/avi"

mkdir -p "$gif_dir" "$avi_dir"

echo "Extracting frames..."

# 1. Decode both into indexed image sequences (TRUE frame index space)
ffmpeg -v error -i "$gif" \
-vf "split=2[base][tmp]; \
[tmp]crop=iw/2:ih:0:0, vflip[left]; \
[base]crop=iw/2:ih:iw/2:0[right]; \
[left][right]hstack=inputs=2" \
"$gif_dir/%06d.png"
ffmpeg -v error -i "$avi" "$avi_dir/%06d.png"

# 2. Determine frame counts
gif_count=$(ls "$gif_dir" | wc -l | tr -d ' ')
avi_count=$(ls "$avi_dir" | wc -l | tr -d ' ')

echo "GIF frames: $gif_count"
echo "AVI frames: $avi_count"

# 3. Choose common length (true index alignment requires this)
N=$(( gif_count < avi_count ? gif_count : avi_count ))

echo "Using N=$N frames for perfect alignment"

# 4. Build trimmed, index-aligned sequences
mkdir -p "$tmpdir/gif_trim" "$tmpdir/avi_trim"

for i in $(seq 1 $N); do
    printf -v idx "%06d" "$i"
    cp "$gif_dir/$idx.png" "$tmpdir/gif_trim/$idx.png"
    cp "$avi_dir/$idx.png" "$tmpdir/avi_trim/$idx.png"
done

# 5. Re-encode with perfect frame pairing
ffmpeg -y \
-framerate "$FPS" -i "$tmpdir/avi_trim/%06d.png" \
-framerate "$FPS" -i "$tmpdir/gif_trim/%06d.png" \
-filter_complex "
[0:v][1:v]scale2ref=-1:ih[avi][gif];
[gif]scale=iw*1.5:ih[gif_stretched];
[avi][gif_stretched]hstack=inputs=2,setpts=PTS*6
" \
-c:v libx264 -pix_fmt yuv420p -crf 18 -movflags +faststart \
"$out"

rm -rf "$tmpdir"

echo "Done: $out"
