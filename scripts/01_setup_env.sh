#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Environment setup for the ctDNA pipeline.
#
# Creates the ctdna2 conda environment from environment.yml if it does not
# already exist, then locates the fgbio JAR and writes scripts/fgbio_env.sh
# for use by 14_align_arms_cd.sh.
#
# fgbio is handled in two ways, in order:
#   1. Conda install (via environment.yml) — preferred; JAR found under
#      $CONDA_PREFIX/share after creation.
#   2. JAR download fallback — used if fgbio was excluded from environment.yml
#      due to gatk4 conflicts or was not found after conda install.
# =============================================================================

SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPTS_DIR/.." && pwd)"
FGBIO_ENV_FILE="$SCRIPTS_DIR/fgbio_env.sh"

FGBIO_VERSION="2.3.0"
FGBIO_JAR_TARGET="$HOME/tools/fgbio/fgbio-${FGBIO_VERSION}.jar"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# System packages — idempotent via apt
# ---------------------------------------------------------------------------
log "Installing system packages..."
sudo apt-get update -qq
sudo apt-get install -y \
  build-essential git curl wget unzip htop tmux \
  python3 python3-pip python3-venv \
  openjdk-17-jdk pigz cloud-guest-utils

# ---------------------------------------------------------------------------
# Miniconda — skip if already installed
# ---------------------------------------------------------------------------
if [[ ! -d "$HOME/miniconda3" ]]; then
  log "Installing Miniconda..."
  wget -q -O "$HOME/Miniconda3-latest-Linux-x86_64.sh" \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash "$HOME/Miniconda3-latest-Linux-x86_64.sh" -b -p "$HOME/miniconda3"
else
  log "Miniconda already installed — skipping."
fi

export PATH="$HOME/miniconda3/bin:$PATH"
eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"

conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r || true

# ---------------------------------------------------------------------------
# ctdna2 conda environment
# ---------------------------------------------------------------------------
ENV_YML="$PROJECT_ROOT/environment.yml"

if conda env list | awk '{print $1}' | grep -qx ctdna2; then
  log "ctdna2 environment already exists — skipping creation."
  log "To recreate: conda env remove -n ctdna2 && bash $0"
else
  if [[ ! -f "$ENV_YML" ]]; then
    echo "ERROR: environment.yml not found at $ENV_YML" >&2
    exit 1
  fi
  log "Creating ctdna2 environment from environment.yml..."
  conda env create -f "$ENV_YML"
  log "ctdna2 environment created."
fi

conda activate ctdna2

# ---------------------------------------------------------------------------
# Locate or download fgbio 2.3.0 JAR
#
# fgbio is NOT in environment.yml — bioconda only packages 1.4.0 and it
# conflicts with gatk4. We always use the upstream release JAR at a fixed
# path. The conda share directory is never searched.
# ---------------------------------------------------------------------------
log "Locating fgbio JAR..."

FGBIO_JAR_PATH="$FGBIO_JAR_TARGET"

if [[ -f "$FGBIO_JAR_PATH" ]]; then
  log "Found fgbio JAR: $FGBIO_JAR_PATH"
else
  log "Downloading fgbio ${FGBIO_VERSION} JAR..."
  mkdir -p "$(dirname "$FGBIO_JAR_PATH")"
  wget -q -O "$FGBIO_JAR_PATH" \
    "https://github.com/fulcrumgenomics/fgbio/releases/download/${FGBIO_VERSION}/fgbio-${FGBIO_VERSION}.jar"
  log "Downloaded: $FGBIO_JAR_PATH"
fi

# Verify
{ java -Xmx512m -jar "$FGBIO_JAR_PATH" --version 2>&1 || true; } | grep -q "Version:" \
  || { echo "ERROR: fgbio JAR did not print expected version string: $FGBIO_JAR_PATH" >&2; exit 1; }
log "fgbio JAR verified."

# Write path for downstream scripts
echo "export FGBIO_JAR=${FGBIO_JAR_PATH}" > "$FGBIO_ENV_FILE"
log "Wrote FGBIO_JAR to $FGBIO_ENV_FILE"

log "Script 01 complete. Activate environment with: conda activate ctdna2"
