SRAID="SRX11133592"
prefetch $SRAID --max-size u
fasterq-dump --split-files --include-technical "$SRAID"sra
