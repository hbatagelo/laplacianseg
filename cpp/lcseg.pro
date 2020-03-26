QT += concurrent
CONFIG += c++17 console
CONFIG -= app_bundle

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

SOURCES += \
        color.cpp \
        commandlineparser.cpp \
        main.cpp \
        seededsegmentationhard.cpp \
        seededsegmentationsoft.cpp \
        seededspsegmentation.cpp \
        slic.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    color.h \
    commandlineparser.h \
    seededsegmentationhard.h \
    seededsegmentationsoft.h \
    seededspsegmentation.h \
    slic.h \
    util.h

RESOURCES +=

win32: {
INCLUDEPATH += "C:\\opencv\\include" # Note: update with the correct OpenCV include path
OPENCV_LIB_PATH = "C:\\opencv\\lib" # Note: update with the correct OpenCV lib path

    debug: {
        LIBS += -L$$OPENCV_LIB_PATH \
                -lopencv_core348d \
                -lopencv_imgproc348d
    }

    release: {
        LIBS += -L$$OPENCV_LIB_PATH \
                -lopencv_core348 \
                -lopencv_imgproc348
    }
}

unix: !macx {
INCLUDEPATH += /usr/include/suitesparse
LIBS += -lcholmod
QT_CONFIG -= no-pkg-config
CONFIG += link_pkgconfig
PKGCONFIG += opencv eigen3
}

unix: macx{
INCLUDEPATH += /usr/local/include
INCLUDEPATH += /usr/local/include/eigen3
LIBS        += -L/usr/local/lib \
               -lopencv_world \
               -lspqr -lumfpack -lcholmod -lmetis -lcolamd -lccolamd -lcamd -lamd -lm -lsuitesparseconfig
}
