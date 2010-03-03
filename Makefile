#----------------------------------------------------------------
#   Parent Makefile for SPLASH
#   This file is just a wrapper for the various sub-makes
#
#   (c) 2007-2010 Daniel Price
#
#----------------------------------------------------------------

.PHONY: splash
splash:
        @cd build; ${MAKE} ${MAKECMDGOALS}

%::
        @cd build; ${MAKE} ${MAKECMDGOALS}

clean:
        @cd build; ${MAKE} clean
