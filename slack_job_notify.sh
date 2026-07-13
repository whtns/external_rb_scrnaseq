#!/bin/bash
# Post a Slack ping with a job's final state. Usage: slack_job_notify.sh <JID> [label]
# Kept as a script file (not an sbatch --wrap) so the JSON quoting can't break.
JID="$1"
LABEL="${2:-job}"
ST=$(sacct -j "$JID" -n -o State | head -1 | xargs)
curl -s -X POST -H 'Content-type: application/json' \
  --data "{\"text\":\"$LABEL ($JID) finished: $ST\"}" "$(cat ~/.slack_webhook)"
