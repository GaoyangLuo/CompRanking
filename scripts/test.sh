#!/bin/bash
set -e
set -m
# default parameters
PREFIX="CompRanking"
THREADS=16
# CONDA_BIN_PATH=~/miniconda/bin
# 获取当前脚本所在目录的上一级目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT_DIR="$(dirname "$SCRIPT_DIR")"

# 读取 YAML 文件中的 conda bin 路径
YAML_FILE="$PARENT_DIR/test_yaml.yaml"

# 使用 Python 解析 YAML
CONDA_BIN_PATH=$(python3 -c "
import yaml
with open('$YAML_FILE', 'r') as f:
    data = yaml.safe_load(f)
print(data['CompRanking']['abs_path_to_conda_bin'])
")

echo $CONDA_BIN_PATH
# 确保获取到路径
if [[ -z "$CONDA_BIN_PATH" ]]; then
    echo "Error: Failed to get Conda bin path from $YAML_FILE"
    exit 1
fi