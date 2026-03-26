#!/bin/bash
# Start Caddy as HTTPS reverse proxy to Django (runserver or gunicorn on :8000).
#
# Prerequisites:
#   1. caddy v2 — https://caddyserver.com/docs/install
#   2. TLS files from ./issue-mkcert.sh (default paths below)
#
# Environment (optional; set in shell or .env):
#   MAS_PUBLIC_HOST   — hostname users open (default: mas.local). Add to /etc/hosts or DNS.
#   MAS_HTTPS_PORT    — listen port (default: 8443). Use 443 only with sudo or cap_net_bind_service.
#   MAS_BACKEND_PORT  — Django port (default: 8000)
#   MAS_TLS_CERT      — path to fullchain/cert PEM
#   MAS_TLS_KEY       — path to private key PEM
#   MAS_BEHIND_HTTPS_PROXY — set to 1 in the environment where Django runs (see below)
#
# Django must be started with proxy-aware settings:
#   export MAS_BEHIND_HTTPS_PROXY=1
#   python manage.py runserver 0.0.0.0:8000
#
# Then open: https://MAS_PUBLIC_HOST:MAS_HTTPS_PORT/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CERT_DIR="${MAS_TLS_CERT_DIR:-$SCRIPT_DIR/certs}"

if [ -f "$PROJECT_ROOT/.env" ]; then
	set -a
	# shellcheck source=/dev/null
	source "$PROJECT_ROOT/.env"
	set +a
fi

export MAS_PUBLIC_HOST="${MAS_PUBLIC_HOST:-mas.local}"
export MAS_HTTPS_PORT="${MAS_HTTPS_PORT:-8443}"
export MAS_BACKEND_HOST="${MAS_BACKEND_HOST:-127.0.0.1}"
export MAS_BACKEND_PORT="${MAS_BACKEND_PORT:-8000}"
export MAS_TLS_CERT="${MAS_TLS_CERT:-$CERT_DIR/mas.pem}"
export MAS_TLS_KEY="${MAS_TLS_KEY:-$CERT_DIR/mas-key.pem}"

if ! command -v caddy &>/dev/null; then
	echo "caddy not found. Install: https://caddyserver.com/docs/install"
	exit 1
fi

if [ ! -f "$MAS_TLS_CERT" ] || [ ! -f "$MAS_TLS_KEY" ]; then
	echo "Missing TLS files: $MAS_TLS_CERT / $MAS_TLS_KEY"
	echo "Run: $SCRIPT_DIR/issue-mkcert.sh"
	exit 1
fi

if [ "${MAS_HTTPS_PORT}" = "443" ] && [ "$(id -u)" -ne 0 ]; then
	echo "Binding to port 443 usually requires sudo. Use MAS_HTTPS_PORT=8443 (default) or run: sudo $0"
	exit 1
fi

echo "HTTPS: https://${MAS_PUBLIC_HOST}:${MAS_HTTPS_PORT}/  →  http://${MAS_BACKEND_HOST}:${MAS_BACKEND_PORT}/"
echo "Ensure Django has MAS_BEHIND_HTTPS_PROXY=1 (and DEVELOPER_MODE=TRUE for dev)."
echo ""

exec caddy run --config "$SCRIPT_DIR/Caddyfile" --adapter caddyfile
