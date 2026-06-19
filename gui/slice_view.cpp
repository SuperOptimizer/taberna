/* slice_view.cpp — see slice_view.h. */
#include "slice_view.h"
#include "volume_source.h"

#include <QGraphicsScene>
#include <QPainter>
#include <QWheelEvent>

SliceView::SliceView(int axis, QWidget *parent) : QGraphicsView(parent), axis_(axis) {
  setScene(new QGraphicsScene(this));
  setBackgroundBrush(Qt::black);
  setDragMode(QGraphicsView::ScrollHandDrag);
  setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
}

void SliceView::setSource(const VolumeSource *src, int lod) {
  src_ = src;
  lod_ = lod;
  if (!src_) return;
  int extent = axis_ == 0 ? src_->nzAt(lod_) : axis_ == 1 ? src_->nyAt(lod_) : src_->nxAt(lod_);
  index_ = extent / 2;
  refresh();
}

void SliceView::setIndex(int index) {
  index_ = index;
  refresh();
}

void SliceView::refresh() {
  if (!src_ || !src_->isOpen()) return;
  src_->readSlice(axis_, index_, lod_, buf_, w_, h_);
  if (w_ <= 0 || h_ <= 0) return;
  image_ = QImage(buf_.data(), w_, h_, w_, QImage::Format_Grayscale8);
  scene()->setSceneRect(0, 0, w_, h_);
  fitInView(scene()->sceneRect(), Qt::KeepAspectRatio);
  viewport()->update();
}

void SliceView::drawBackground(QPainter *p, const QRectF & /*rect*/) {
  if (!image_.isNull()) p->drawImage(QPointF(0, 0), image_);
}

void SliceView::wheelEvent(QWheelEvent *e) {
  if (e->modifiers() & Qt::ControlModifier) {  // ctrl+wheel = zoom
    double s = e->angleDelta().y() > 0 ? 1.15 : 1.0 / 1.15;
    scale(s, s);
    return;
  }
  if (!src_) return;
  int extent = axis_ == 0 ? src_->nzAt(lod_) : axis_ == 1 ? src_->nyAt(lod_) : src_->nxAt(lod_);
  int step = e->angleDelta().y() > 0 ? 1 : -1;
  index_ += step;
  if (index_ < 0) index_ = 0;
  if (index_ >= extent) index_ = extent - 1;
  refresh();
}

void SliceView::resizeEvent(QResizeEvent *e) {
  QGraphicsView::resizeEvent(e);
  if (!image_.isNull()) fitInView(scene()->sceneRect(), Qt::KeepAspectRatio);
}
