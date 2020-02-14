SOURCES += \
    mainwindowCD.cpp \
    qcustomplotCD.cpp \
    IS_cell.cpp

HEADERS += \
    mainwindowCD.h \
    qcustomplotCD.h \
    IS_cell.h

QMAKE_CXXFLAGS += -std=c++11 -O3

QT += printsupport

TARGET = CD2_corolla
FORMS += \
    mainwindowCD.ui

