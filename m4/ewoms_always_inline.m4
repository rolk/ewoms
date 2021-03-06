# Check whether the function attibute 'always_inline' is available and sets
# the EWOMS_ALWAYS_INLINE macro accordingly
AC_DEFUN([EWOMS_CHECK_ALWAYS_INLINE],[
  AC_MSG_CHECKING([for __attribute__((always_inline))])
  AC_LANG_PUSH([C++])
  AC_TRY_COMPILE([void fn() __attribute__((always_inline)); void fn() {} ],[],
                 [EWOMS_ALWAYS_INLINE="__attribute__((always_inline))"
                  AC_MSG_RESULT(yes)],
                 [EWOMS_ALWAYS_INLINE=""
                  AC_MSG_RESULT(no)])
  AC_LANG_POP([C++])
  AC_DEFINE_UNQUOTED(EWOMS_ALWAYS_INLINE, $EWOMS_ALWAYS_INLINE,
                     [USE WITH CARE: Forces a function to be inlined even for non-optimized builds])
])
