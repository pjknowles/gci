#-------------------------------------------------
#
# Project created by QtCreator 2012-03-06T06:11:54
#
#-------------------------------------------------

CONFIG -= qt
QT       -= core

QT       -= gui

TARGET = gci.exe
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

#INCLUDEPATH += /opt/intel/composerxe/mkl/include
#LIBS += -L/opt/intel/composerxe/mkl/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#LIBS += -L/usr/local/Cellar/google-perftools/2.1/lib -lprofiler
QMAKE_CXXFLAGS_X86_64 -= -Xarch_x86_64
QMAKE_LFLAGS_X86_64 -= -Xarch_x86_64
QMAKE_CFLAGS_X86_64 -= -Xarch_x86_64
QMAKE_LINK = $$QMAKE_CXX # why would you want it different?

contains(QMAKE_CXX,mpic++) {
message(parallel build)
#QMAKE_CXXFLAGS=
GA=/usr/local
PPIDD=/usr/local
INCLUDEPATH+=$$PPIDD/src $$PPIDD/include
DEFINES += GCI_PARALLEL
LIBS += -L$$PPIDD -lppidd -L$$GA/lib -lga -larmci -lblas
}
INCLUDEPATH+=FCIdump

SOURCES += gci.cpp \
    gciHamiltonian.cpp \
    gciState.cpp \
    #gciNode.cpp \
    #gciExcitation.cpp \
    gciDeterminant.cpp \
    gciString.cpp \
    gciWavefunction.cpp \
    FCIdump/FCIdump.cpp \
    gciStringSet.cpp \
    gciExcitationSet.cpp \
    gciOperator.cpp \
    gciPrintable.cpp \
    gciSymmetrySpace.cpp \
    gciOrbitalSpace.cpp \
    gciTransitionDensity.cpp \
    gciFile.cpp \
    gciMolpro.cpp \
    gciRun.cpp \
    Profiler.cpp

HEADERS += gci.h \
    gciHamiltonian.h \
    gciState.h \
    #gciNode.h \
    #gciExcitation.h \
    gciDeterminant.h \
    gciString.h \
    gciWavefunction.h \
    FCIdump/FCIdump.h \
    gciStringSet.h \
    gciExcitationSet.h \
    gciOperator.h \
    gciPrintable.h \
    gciSymmetrySpace.h \
    gciOrbitalSpace.h \
    gciTransitionDensity.h \
    gciFile.h \
    gciMolpro.h \
    gciRun.h \
    Profiler.h ProfilerC.h

html.target = $$PWD/html
html.commands = (cd $$PWD ; doxygen $$PWD/Doxyfile;)
html.depends = $(OBJECTS)
dox.target = html
dox.depends = $$PWD/html
QMAKE_EXTRA_TARGETS += dox html all

all.target = all
all.depends += gci.fcidump
all.depends += html

input.target = gci.fcidump
input.commands = ln -s $$PWD/gci.fcidump gci.fcidump
input.depends =
QMAKE_EXTRA_TARGETS += input






