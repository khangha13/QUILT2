#!/bin/bash
# QUILT2 SLURM resource defaults (mirrors Step1C pattern).
# Edit this file to change resource requests; environment.sh is for paths,
# tools, and behavior toggles.

# Account/partition
QUILT2_ACCOUNT="a_qaafi_cas"
QUILT2_PARTITION="general"
QUILT2_QOS=""

# Resources
QUILT2_NODES="1"
QUILT2_NTASKS="1"
QUILT2_CPUS_PER_TASK="2"
QUILT2_PHASE2_CPUS_PER_TASK="2"
QUILT2_MEMORY="12G"
QUILT2_TIME_LIMIT="72:00:00"
QUILT2_MASTER_TIME_LIMIT="336:00:00"

# Array limit (max tasks submitted even if manifest longer)
QUILT2_ARRAY_MAX="0" # 0 = no cap
# Default to epyc4 nodes to avoid bcftools ISA issues on older hardware
QUILT2_CONSTRAINT="epyc4"

get_quilt2_config() {
    cat <<EOF
ACCOUNT=${QUILT2_ACCOUNT}
PARTITION=${QUILT2_PARTITION}
QOS=${QUILT2_QOS}
NODES=${QUILT2_NODES}
NTASKS=${QUILT2_NTASKS}
CPUS=${QUILT2_CPUS_PER_TASK}
PHASE2_CPUS=${QUILT2_PHASE2_CPUS_PER_TASK}
MEMORY=${QUILT2_MEMORY}
TIME=${QUILT2_TIME_LIMIT}
MASTER_TIME=${QUILT2_MASTER_TIME_LIMIT}
ARRAY_MAX=${QUILT2_ARRAY_MAX}
CONSTRAINT=${QUILT2_CONSTRAINT}
EOF
}
