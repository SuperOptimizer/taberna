/* slice_view.h — one orthogonal MPR slice pane, mirroring VC3D's CVolumeViewerView:
 * a QGraphicsView that blits a CPU-rendered slice QImage in drawBackground (so the
 * scene graph is free for overlays later) and scrolls through slices on the wheel.
 * Data comes from VolumeSource (axis 0=axial/1=coronal/2=sagittal).
 */
#pragma once
#include <QGraphicsView>
#include <QImage>
#include <cstdint>
#include <vector>

class VolumeSource;

class SliceView : public QGraphicsView {
  Q_OBJECT
public:
  explicit SliceView(int axis, QWidget *parent = nullptr);
  void setSource(const VolumeSource *src, int lod);
  void setIndex(int index);
  int  index() const { return index_; }

protected:
  void drawBackground(QPainter *p, const QRectF &rect) override;
  void wheelEvent(QWheelEvent *e) override;
  void resizeEvent(QResizeEvent *e) override;

private:
  void refresh();

  const VolumeSource  *src_ = nullptr;
  int                  axis_ = 0;
  int                  lod_ = 0;
  int                  index_ = 0;
  int                  w_ = 0, h_ = 0;
  std::vector<uint8_t> buf_;   // backing store for image_
  QImage               image_;
};
