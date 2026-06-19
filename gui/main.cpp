/* main.cpp — taberna Qt+VTK 2x2 volume viewer entry point. */
#include "main_window.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>

int main(int argc, char **argv) {
  // must be set before any QVTKOpenGLNativeWidget is created
  QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

  QApplication app(argc, argv);
  MainWindow w;
  w.show();
  if (argc > 1) w.openArchive(QString::fromLocal8Bit(argv[1]));
  return app.exec();
}
