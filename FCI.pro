#-------------------------------------------------
#
# Project created by QtCreator 2012-03-06T06:11:54
#
#-------------------------------------------------

QT       -= core

QT       -= gui

TARGET = FCI
CONFIG   += console
CONFIG   -= app_bundle

#INCLUDEPATH += /usr/local/boost

TEMPLATE = app


SOURCES += main.cpp \
    FCIHamiltonian.cpp \
    FCIParameters.cpp \
    FCINode.cpp \
    FCIExcitation.cpp \
    FCIDeterminant.cpp \
    FCIString.cpp

HEADERS += \
    FCIHamiltonian.h \
    FCIParameters.h \
    FCINode.h \
    FCIExcitation.h \
    FCIDeterminant.h \
    FCIString.h





