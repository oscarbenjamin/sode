#!/usr/bin/env python

from sode.cysode import CYSODE
from sode.script import Script

script = Script(CYSODE)

if __name__ == "__main__":
    import sys
    script.main(argv=sys.argv[1:])
