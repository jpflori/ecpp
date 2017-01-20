# Libtool -version-info CURRENT[:REVISION[:AGE]] for libmpfr.la
#
# 1. No interfaces changed, only implementations (good):
#    ==> Increment REVISION.
# 2. Interfaces added, none removed (good):
#    ==> Increment CURRENT, increment AGE, set REVISION to 0.
# 3. Interfaces removed or changed (BAD, breaks upward compatibility):
#    ==> Increment CURRENT, set AGE and REVISION to 0.

LIBECPP_LT_CURRENT    = 0
LIBECPP_LT_REVISION   = 0
LIBECPP_LT_AGE        = 0

include_HEADERS = ecpp.h

lib_LTLIBRARIES = libecpp.la

libecpp_la_SOURCES = primes.c indices.c disc.c sqrt_mod.c cornacchia2.c factor.c step.c

libecpp_la_LDFLAGS = $(ECPP_LDFLAGS) $(LIBECPP_LDFLAGS) \
  -version-info $(LIBECPP_LT_CURRENT):$(LIBECPP_LT_REVISION):$(LIBECPP_LT_AGE)
