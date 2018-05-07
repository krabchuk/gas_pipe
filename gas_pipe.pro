TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CFLAGS += -Wall -Wextra -Werror

SOURCES += main.c \
    laspack/eigenval.c \
    laspack/errhandl.c \
    laspack/factor.c \
    laspack/itersolv.c \
    laspack/matrix.c \
    laspack/mlsolv.c \
    laspack/operats.c \
    laspack/precond.c \
    laspack/qmatrix.c \
    laspack/rtc.c \
    laspack/vector.c \
    laspack/eigenval.c \
    laspack/errhandl.c \
    laspack/factor.c \
    laspack/itersolv.c \
    laspack/matrix.c \
    laspack/mlsolv.c \
    laspack/operats.c \
    laspack/precond.c \
    laspack/qmatrix.c \
    laspack/rtc.c \
    laspack/vector.c \
    main.c \
    solver.c \
    parser.c \
    neighbors_and_type_init.c

HEADERS += \
    laspack/copyrght.h \
    laspack/eigenval.h \
    laspack/elcmp.h \
    laspack/errhandl.h \
    laspack/factor.h \
    laspack/itersolv.h \
    laspack/lastypes.h \
    laspack/matrix.h \
    laspack/mlsolv.h \
    laspack/operats.h \
    laspack/precond.h \
    laspack/qmatrix.h \
    laspack/rtc.h \
    laspack/vector.h \
    laspack/version.h \
    laspack_include.h \
    laspack/copyrght.h \
    laspack/eigenval.h \
    laspack/elcmp.h \
    laspack/errhandl.h \
    laspack/factor.h \
    laspack/itersolv.h \
    laspack/lastypes.h \
    laspack/matrix.h \
    laspack/mlsolv.h \
    laspack/operats.h \
    laspack/precond.h \
    laspack/qmatrix.h \
    laspack/rtc.h \
    laspack/vector.h \
    laspack/version.h \
    laspack_include.h \
    solver.h \
    parser.h \
    neighbors_and_type_init.h
