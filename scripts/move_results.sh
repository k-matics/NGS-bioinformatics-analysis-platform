#!/usr/bin/env bash
# Move the platform's run outputs to a new results folder, safely.
#
#   scripts/move_results.sh <old_results_dir> <new_results_dir> [--delete-old]
#
# Steps (stops at the first problem; nothing is deleted unless every step passed):
#   1. refuses if a workflow is running/queued, or the new folder is not empty
#   2. copies with rsync -a (permissions, times and symlinks kept)
#   3. verifies EVERY file: size + SHA-256 of source vs copy; writes a manifest into the new folder
#   4. re-points absolute symlinks that pointed inside the old folder
#   5. backs up the database, then rewrites stored paths old -> new (sample files, run folders,
#      step outputs / logs / summaries). Command history and event messages are left unchanged:
#      they record what actually ran.
#   6. writes NGS_RESULTS_DIR into config/local.env
#   7. only with --delete-old: removes the old folder
# Services must be stopped before running this (scripts/stop_services.sh).
set -euo pipefail

OLD="$(realpath -e "${1:?old results dir}")"
NEW="${2:?new results dir}"
DELETE="${3:-}"
APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="$APP_DIR/database/ngs_platform.db"
STAMP="$(date +%Y%m%d-%H%M%S)"

[[ -d "$NEW" && -w "$NEW" ]] || { echo "ERROR: $NEW does not exist or is not writable"; exit 1; }
NEW="$(realpath -e "$NEW")"
[[ "$OLD" != "$NEW" ]] || { echo "ERROR: same folder"; exit 1; }
if [[ -n "$(ls -A "$NEW")" ]]; then echo "ERROR: $NEW is not empty — refusing to overwrite"; exit 1; fi
for pid in instance/worker.pid instance/web.pid; do
  if [[ -f "$APP_DIR/$pid" ]] && kill -0 "$(cat "$APP_DIR/$pid")" 2>/dev/null; then
    echo "ERROR: services are running ($pid) — run scripts/stop_services.sh first"; exit 1
  fi
done
ACTIVE=$(sqlite3 "$DB" "select count(*) from workflow_runs where status in ('running','queued','pending','paused')")
[[ "$ACTIVE" == 0 ]] || { echo "ERROR: $ACTIVE workflow(s) not finished"; exit 1; }

need=$(du -sb "$OLD" | cut -f1); free=$(df -B1 --output=avail "$NEW" | tail -1)
echo "Need $(numfmt --to=iec "$need"), free on target $(numfmt --to=iec "$free")"
MARGIN_GB="${MOVE_MARGIN_GB:-10}"
(( free > need + MARGIN_GB*1024**3 )) || { echo "ERROR: not enough space on target (keeps ${MARGIN_GB} GB margin)"; exit 1; }

echo "[1/6] Copying $OLD -> $NEW"
rsync -a --info=stats1 "$OLD"/ "$NEW"/

echo "[2/6] Verifying every file (size + SHA-256)"
MAN="$NEW/MIGRATION_MANIFEST_$STAMP.tsv"
printf "relative_path\tsize_bytes\tsha256\n" > "$MAN"
bad=0; n=0
while IFS= read -r -d '' f; do
  rel="${f#"$OLD"/}"
  s1=$(stat -c %s "$f"); s2=$(stat -c %s "$NEW/$rel" 2>/dev/null || echo missing)
  h1=$(sha256sum "$f" | cut -d' ' -f1); h2=$(sha256sum "$NEW/$rel" 2>/dev/null | cut -d' ' -f1 || true)
  if [[ "$s1" != "$s2" || "$h1" != "$h2" ]]; then echo "  MISMATCH: $rel"; bad=$((bad+1)); fi
  printf "%s\t%s\t%s\n" "$rel" "$s1" "$h1" >> "$MAN"; n=$((n+1))
done < <(find "$OLD" -type f -print0)
nl_old=$(find "$OLD" -type l | wc -l); nl_new=$(find "$NEW" -type l | wc -l)
echo "  $n files checked, $bad mismatch(es); symlinks $nl_old -> $nl_new"
[[ "$bad" == 0 && "$nl_old" == "$nl_new" ]] || { echo "ERROR: verification failed — old folder untouched"; exit 1; }

echo "[3/6] Re-pointing symlinks that pointed inside the old folder"
while IFS= read -r -d '' l; do
  t="$(readlink "$l")"
  if [[ "$t" == "$OLD"/* ]]; then ln -sfn "$NEW/${t#"$OLD"/}" "$l"; echo "  ${l#"$NEW"/}"; fi
done < <(find "$NEW" -type l -print0)

echo "[4/6] Updating stored paths in the database (backup first)"
cp -p "$DB" "$DB.before-results-move-$STAMP"
sqlite3 "$DB" <<SQL
BEGIN;
UPDATE sample_files   SET path       = replace(path,       '$OLD/', '$NEW/') WHERE path       LIKE '$OLD/%';
UPDATE workflow_runs  SET output_dir = replace(output_dir, '$OLD/', '$NEW/') WHERE output_dir LIKE '$OLD/%';
UPDATE step_runs      SET log_path   = replace(log_path,   '$OLD/', '$NEW/') WHERE log_path   LIKE '$OLD/%';
UPDATE step_runs      SET outputs    = replace(outputs,    '$OLD/', '$NEW/') WHERE outputs    LIKE '%$OLD/%';
UPDATE step_runs      SET summary    = replace(summary,    '$OLD/', '$NEW/') WHERE summary    LIKE '%$OLD/%';
COMMIT;
SQL
left=$(sqlite3 "$DB" "select (select count(*) from sample_files where path like '$OLD/%') + (select count(*) from workflow_runs where output_dir like '$OLD/%') + (select count(*) from step_runs where log_path like '$OLD/%' or outputs like '%$OLD/%' or summary like '%$OLD/%')")
[[ "$left" == 0 ]] || { echo "ERROR: $left stored path(s) still point to the old folder"; exit 1; }
echo "  database backup: $DB.before-results-move-$STAMP"

echo "[5/6] Writing NGS_RESULTS_DIR to config/local.env"
ENVF="$APP_DIR/config/local.env"
touch "$ENVF"
grep -v '^NGS_RESULTS_DIR=' "$ENVF" > "$ENVF.tmp" || true
echo "NGS_RESULTS_DIR=$NEW" >> "$ENVF.tmp"
mv "$ENVF.tmp" "$ENVF"
cat "$ENVF"

if [[ "$DELETE" == "--delete-old" ]]; then
  echo "[6/6] Removing the old folder (all files verified)"
  rm -rf -- "$OLD"
else
  echo "[6/6] Old folder kept: $OLD (re-run is not needed; delete it yourself when satisfied)"
fi
echo "Done. Start the services: scripts/start_services.sh"
