#!/bin/bash
# Generate a trusted (mkcert CA) certificate for MAS HTTPS reverse proxy.
# Install mkcert first: https://github.com/FiloSottile/mkcert
#   Linux: sudo apt install mkcert libnss3-tools  OR  brew install mkcert
# Then run once: mkcert -install   (installs the local CA on this machine)
#
# Usage:
#   ./issue-mkcert.sh                    # cert for MAS_PUBLIC_HOST (default mas.local) only
#   ./issue-mkcert.sh 192.168.1.42       # include LAN IP in SAN (recommended for IP access)
#   MAS_PUBLIC_HOST=mas.lab ./issue-mkcert.sh 10.0.0.5

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CERT_DIR="${MAS_TLS_CERT_DIR:-$SCRIPT_DIR/certs}"
HOST="${MAS_PUBLIC_HOST:-mas.local}"

if ! command -v mkcert &>/dev/null; then
	echo "mkcert not found. Install it: https://github.com/FiloSottile/mkcert"
	exit 1
fi

mkdir -p "$CERT_DIR"
CERT_FILE="$CERT_DIR/mas.pem"
KEY_FILE="$CERT_DIR/mas-key.pem"

echo "Writing certificate for: $HOST $*"
echo "Output: $CERT_FILE / $KEY_FILE"

mkcert -cert-file "$CERT_FILE" -key-file "$KEY_FILE" "$HOST" "$@"

echo "Done. Install the mkcert CA on other machines (run 'mkcert -CAROOT' to find the root CA) so browsers trust this cert."
