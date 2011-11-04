#!/usr/bin/env python

import sys

from sode import Script, Weiner

if __name__ == "__main__":
    s = Script(Weiner)
    s.main(sys.argv[1:])
