#!/usr/bin/env bash
# Build the solver and launch the NACA 0012, Re = 1000, PIMPLE cases in the background.
# Usage (from code/):
#   ./run_all.sh sweep        # 8 angles of attack (0, 4, 8, 10, 11, 12, 15, 20 deg), about 2.5 h with all in parallel
#   ./run_all.sh grid         # medium grid 525x375 at 10 deg (coarse level = sweep aoa10), about 3-4 h
#   ./run_all.sh video        # 6 animation runs, a VTK frame every 0.1 for t >= 40, about 2 h
#   ./run_all.sh all          # everything above
# Runs are detached with nohup, so they survive closing the terminal (not sleep/shutdown of the PC).
# Logs: runs/logs/<case>.log. Then post-process with:
#   python3 scripts/postprocess_report.py runs ../results
#   scripts/render_videos.sh
set -euo pipefail
cd "$(dirname "$0")"
what=${1:-sweep}
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build build -j"$(nproc)" >/dev/null
mkdir -p runs/logs
launch() { for cfg in "$@"; do
  name=$(basename "$cfg" .cfg)
  setsid nohup stdbuf -oL -eL ./build/airfoil "$cfg" > "runs/logs/$name.log" 2>&1 < /dev/null &
  echo "started $name"
done; }
case $what in
  sweep) launch cases/sweep/*.cfg ;;
  grid)  launch cases/grid_study/*.cfg ;;
  video) launch cases/video/*.cfg ;;
  all)   launch cases/sweep/*.cfg cases/grid_study/*.cfg cases/video/*.cfg ;;
  *) echo "unknown option: $what (sweep | grid | video | all)"; exit 1 ;;
esac
echo "progress: tail -n1 runs/logs/*.log"
