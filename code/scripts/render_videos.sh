#!/usr/bin/env bash
# Render all vortex-shedding videos from the animation runs (cases/video/*.cfg -> runs/video/).
# Usage (from code/): scripts/render_videos.sh [python]   (python needs numpy, matplotlib >= 3.10 for the 'berlin' map, imageio-ffmpeg or ffmpeg)
set -euo pipefail
PY=${1:-python3}
cd "$(dirname "$0")/.."
OUT=../results/videos; mkdir -p "$OUT"
R=runs/video
for a in 8 10 11 12 15 20; do
  "$PY" scripts/make_vortex_video.py $R/naca0012_Re1000_aoa$a -o $OUT/naca0012_Re1000_aoa$a.mp4 &
done
"$PY" scripts/make_vortex_video.py $R/naca0012_Re1000_aoa{8,10,11,12,15,20} --montage -o $OUT/naca0012_Re1000_all_angles.mp4 &
"$PY" scripts/make_vortex_video.py $R/naca0012_Re1000_aoa{10,12,15,20} --montage -o $OUT/naca0012_Re1000_shedding_4panel.mp4 &
wait
ls -lh "$OUT"
