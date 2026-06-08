#!/usr/bin/env bash
# Run the 01-06 definition smoke tests via wolframscript.
#
#   ./run_tests.sh            # both detectors (TES, MKID)
#   ./run_tests.sh TES        # one detector
#
# Reports land in  Mathematica/output/<DET>/test/report.txt

set -u
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# locate wolframscript (PATH, then the macOS app bundle)
ws="$(command -v wolframscript 2>/dev/null || true)"
if [[ -z "$ws" ]]; then
  for c in /Applications/Wolfram.app/Contents/MacOS/wolframscript \
           /Applications/Mathematica.app/Contents/MacOS/wolframscript; do
    [[ -x "$c" ]] && ws="$c" && break
  done
fi
if [[ -z "$ws" ]]; then
  echo "error: wolframscript not found (PATH or /Applications/*.app)" >&2
  exit 127
fi

dets=("$@")
[[ ${#dets[@]} -eq 0 ]] && dets=(TES MKID)

rc=0
for det in "${dets[@]}"; do
  echo "======================================================================"
  "$ws" -file "$here/test_definitions.wls" "$det" || rc=1
  echo
done

exit "$rc"
