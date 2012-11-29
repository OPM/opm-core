# - Turn on optimizations based on build type

include (AddOptions)

# if we are building a debug target, then disable all optimizations
# otherwise, turn them on. indicate to the code what we have done
# so it can turn on assertions etc.
if (CMAKE_COMPILER_IS_GNUCXX)
  add_options (
	ALL_LANGUAGES
	"Debug"
	"-O0" "-DDEBUG"
	)
  add_options (
	ALL_LANGUAGES
	"Release;RelWithDebInfo;MinSizeRel"
	"-O3" "-DNDEBUG" "-flto"
	)
endif (CMAKE_COMPILER_IS_GNUCXX)
