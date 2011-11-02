#bash

SODEDIR="$(pwd)"
if [ -n "${PYTHONPATH}" ]; then
    PYTHONPATH="${PYTHONPATH}:${SODEDIR}"
else
    PYTHONPATH="${SODEDIR}"
fi

export PYTHONPATH
