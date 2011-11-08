# Copyright Oscar Benjamin 2011
#
# This script is run when sode is called at the cmd line

from sode import Weiner, Script

weinerscript = Script(Weiner)

if __name__ == "__main__":
    import sys
    weinerscript.main(sys.argv[1:])
