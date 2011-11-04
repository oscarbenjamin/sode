#bash

SODEDIR="$(pwd)"
if [ -n "${PYTHONPATH}" ]; then
    PYTHONPATH="${PYTHONPATH}:${SODEDIR}"
else
    PYTHONPATH="${SODEDIR}"
fi

PATH="${PATH}:${SODEDIR}/bin"

export PYTHONPATH
