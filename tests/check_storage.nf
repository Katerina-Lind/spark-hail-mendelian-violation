nextflow.enable.dsl = 2

process CHECK_NVME {
    tag 'nvme-storage-test'

    cpus 2
    memory '4 GB'
    container 'ubuntu:24.04'
    scratch false

    script:
    '''
    set -eu

    echo "===== HOST ====="
    hostname

    echo
    echo "===== ENVIRONMENT ====="
    env | sort | grep -E \
      'TMPDIR|TMP|SPARK|FUSION|NXF' || true

    echo
    echo "===== FILESYSTEMS ====="
    df -hT

    echo
    echo "===== TMPDIR ====="
    echo "TMPDIR=${TMPDIR:-UNSET}"

    LOCAL_TMP="${TMPDIR:-/tmp}"
    echo "Local temporary base: $LOCAL_TMP"

    echo
    echo "===== MOUNT DETAILS ====="
    df -PT "$LOCAL_TMP"
    df -PT /tmp
    cat /proc/self/mountinfo

    echo
    echo "===== WRITE TEST ====="
    mkdir -p "$LOCAL_TMP"
    TESTDIR=$(mktemp -d "$LOCAL_TMP/nvme-check.XXXXXX")
    trap 'rm -rf "$TESTDIR"' EXIT

    dd if=/dev/zero \
       of="$TESTDIR/test.bin" \
       bs=1M count=64 conv=fsync

    echo
    echo "===== TEST FILE LOCATION ====="
    df -hT "$TESTDIR"
    ls -lh "$TESTDIR/test.bin"

    echo
    echo "===== STORAGE TEST COMPLETE ====="
    '''
}

workflow {
    CHECK_NVME()
}
