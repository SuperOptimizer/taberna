/* main_window.h — the 2x2 viewer, mirroring VC3D's QMdiArea-of-4-subwindows
 * layout: three Qt-native ortho slice panes (axial / coronal / sagittal) plus the
 * new VTK volume-render pane. Data flows from a single VolumeSource (matter-
 * compressor) into all four. Trimmed hard from VC3D — no docks, no overlays yet.
 */
#pragma once
#include <QMainWindow>
#include "volume_source.h"

class SliceView;
class VolumeView;

class MainWindow : public QMainWindow {
  Q_OBJECT
public:
  MainWindow();
  void openArchive(const QString &path);

private slots:
  void onOpen();

private:
  int sliceLod() const;
  int volumeLod() const;

  VolumeSource source_;
  SliceView   *axial_ = nullptr;
  SliceView   *coronal_ = nullptr;
  SliceView   *sagittal_ = nullptr;
  VolumeView  *volume_ = nullptr;
};
