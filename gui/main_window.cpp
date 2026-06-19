/* main_window.cpp — see main_window.h. */
#include "main_window.h"
#include "slice_view.h"
#include "volume_view.h"

#include <QFileDialog>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QMenuBar>
#include <QMessageBox>
#include <QStatusBar>
#include <cstdio>
#include <vector>

static QMdiSubWindow *addPane(QMdiArea *mdi, QWidget *w, const QString &title) {
  QMdiSubWindow *sub = mdi->addSubWindow(w);
  sub->setWindowTitle(title);
  sub->setAttribute(Qt::WA_DeleteOnClose, false);
  return sub;
}

MainWindow::MainWindow() {
  auto *mdi = new QMdiArea(this);
  setCentralWidget(mdi);

  axial_ = new SliceView(0);
  coronal_ = new SliceView(1);
  sagittal_ = new SliceView(2);
  volume_ = new VolumeView();

  addPane(mdi, axial_, tr("Axial (XY)"));
  addPane(mdi, coronal_, tr("Coronal (XZ)"));
  addPane(mdi, sagittal_, tr("Sagittal (YZ)"));
  addPane(mdi, volume_, tr("Volume (VTK)"));
  mdi->tileSubWindows();

  auto *file = menuBar()->addMenu(tr("&File"));
  file->addAction(tr("&Open .mca…"), this, &MainWindow::onOpen);
  file->addAction(tr("&Quit"), this, &QWidget::close);

  auto *help = menuBar()->addMenu(tr("&Help"));
  help->addAction(tr("&GPU info"), this, [this] {
    QMessageBox::information(this, tr("GPU info"), volume_->gpuInfo());
  });

  setWindowTitle(tr("taberna viewer"));
  resize(1200, 900);
  statusBar()->showMessage(tr("Open an .mca archive (File ▸ Open)"));
}

void MainWindow::onOpen() {
  QString path = QFileDialog::getOpenFileName(this, tr("Open matter-compressor archive"),
                                              QString(), tr("mca archives (*.mca);;all (*)"));
  if (!path.isEmpty()) openArchive(path);
}

int MainWindow::sliceLod() const {
  int maxd = std::max({source_.nx(), source_.ny(), source_.nz()});
  int lod = 0;
  while ((maxd >> lod) > 2048 && lod < source_.nlods() - 1) lod++;
  return lod;
}

int MainWindow::volumeLod() const {
  // pick a LOD whose whole volume fits comfortably for GPU upload (~256^3)
  int lod = 0;
  while (lod < source_.nlods() - 1) {
    long vox = (long)source_.nxAt(lod) * source_.nyAt(lod) * source_.nzAt(lod);
    if (vox <= 256L * 256 * 256) break;
    lod++;
  }
  return lod;
}

void MainWindow::openArchive(const QString &path) {
  if (!source_.open(path.toStdString())) {
    QMessageBox::warning(this, tr("Open failed"), tr("Could not open %1").arg(path));
    return;
  }
  int slod = sliceLod(), vlod = volumeLod();
  axial_->setSource(&source_, slod);
  coronal_->setSource(&source_, slod);
  sagittal_->setSource(&source_, slod);

  int dz = source_.nzAt(vlod), dy = source_.nyAt(vlod), dx = source_.nxAt(vlod);
  std::vector<uint8_t> vol((size_t)dz * dy * dx, 0);
  source_.readRegion(vlod, 0, 0, 0, dz, dy, dx, vol.data());
  volume_->setVolume(vol.data(), dz, dy, dx);

  setWindowTitle(tr("taberna viewer — %1  [%2x%3x%4, %5 LODs, q=%6]")
                     .arg(path).arg(source_.nx()).arg(source_.ny()).arg(source_.nz())
                     .arg(source_.nlods()).arg(source_.quality()));
  statusBar()->showMessage(tr("slice LOD %1, volume LOD %2").arg(slod).arg(vlod));
  std::fprintf(stderr, "%s\n", volume_->gpuInfo().toLocal8Bit().constData());
}
