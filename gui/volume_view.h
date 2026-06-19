/* volume_view.h — the NEW pane VC3D does not have: a GPU volume renderer.
 * A QVTKOpenGLNativeWidget driving a vtkSmartVolumeMapper over a vtkImageData we
 * fill from a VolumeSource subvolume. setVolume() uploads a dense u8 (dz,dy,dx)
 * block and renders it with a simple grayscale transfer function.
 */
#pragma once
#include <QString>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkColorTransferFunction.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPiecewiseFunction.h>
#include <vtkRenderer.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkVolume.h>
#include <vtkVolumeProperty.h>
#include <cstdint>

class VolumeView : public QVTKOpenGLNativeWidget {
  Q_OBJECT
public:
  explicit VolumeView(QWidget *parent = nullptr);
  // Upload a dense u8 region (z-major, x-fastest) and render it.
  void setVolume(const uint8_t *data, int dz, int dy, int dx);

  // GPU capability report: renderer, GL version, 3D-texture cap, VRAM (if the
  // driver exposes it), and the resulting max u8 cube edge we can render. This is
  // what determines the largest volume the current GPU can hold (see the texture
  // cap vs. memory discussion). Requires a live GL context (call after show()).
  QString gpuInfo();

  // --- adjustable render controls (driven by the VolumeControls dock) ---
  void setWindowLevel(double window, double level);  // CT-style contrast
  void setColormap(int preset);   // 0 gray, 1 bone, 2 hot, 3 viridis
  void setGradientOpacity(bool on);
  void setShade(bool on);
  void setBlendMode(int mode);    // 0 composite,1 MIP,2 MinIP,3 average,4 isosurface
  void setIsoValue(double v);
  void setCrop(bool enable, double fx0, double fx1, double fy0, double fy1,
               double fz0, double fz1);  // fractions in [0,1]
  int dimX() const { return dx_; }
  int dimY() const { return dy_; }
  int dimZ() const { return dz_; }

private:
  void rebuildColorOpacity();
  vtkNew<vtkGenericOpenGLRenderWindow> renderWindow_;
  vtkNew<vtkRenderer>                  renderer_;
  vtkNew<vtkSmartVolumeMapper>         mapper_;
  vtkNew<vtkVolume>                    volume_;
  vtkNew<vtkImageData>                 image_;
  vtkNew<vtkColorTransferFunction>     color_;
  vtkNew<vtkPiecewiseFunction>         opacity_;
  vtkNew<vtkVolumeProperty>            property_;
  vtkNew<vtkPiecewiseFunction>         gradOpacity_;
  bool                                 added_ = false;  // volume actor in renderer yet?
  double level_ = 128.0, window_ = 255.0, iso_ = 128.0;
  int    colormap_ = 0, blend_ = 0;
  int    dx_ = 0, dy_ = 0, dz_ = 0;
};
