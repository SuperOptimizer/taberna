/* volume_view.cpp — see volume_view.h. */
#include "volume_view.h"

#include <QByteArray>
#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vtkGPUVolumeRayCastMapper.h>

// GL enums we query via QOpenGLFunctions (guard if Qt's GL headers omit them).
#ifndef GL_MAX_3D_TEXTURE_SIZE
#define GL_MAX_3D_TEXTURE_SIZE 0x8073
#endif
#ifndef GL_MAX_TEXTURE_SIZE
#define GL_MAX_TEXTURE_SIZE 0x0D33
#endif
#ifndef GL_VENDOR
#define GL_VENDOR 0x1F00
#endif
#ifndef GL_RENDERER
#define GL_RENDERER 0x1F01
#endif
#ifndef GL_VERSION
#define GL_VERSION 0x1F02
#endif
// NVX_gpu_memory_info enums (NVIDIA) — values in KiB.
#ifndef GL_GPU_MEMORY_INFO_TOTAL_AVAILABLE_MEMORY_NVX
#define GL_GPU_MEMORY_INFO_TOTAL_AVAILABLE_MEMORY_NVX 0x9048
#define GL_GPU_MEMORY_INFO_CURRENT_AVAILABLE_VIDMEM_NVX 0x9049
#endif
// ATI_meminfo (AMD) — free memory in KiB (param 0 = total pool).
#ifndef GL_TEXTURE_FREE_MEMORY_ATI
#define GL_TEXTURE_FREE_MEMORY_ATI 0x87FC
#endif

static long sys_mem_available_kb() {
  FILE *f = std::fopen("/proc/meminfo", "r");
  if (!f) return -1;
  char key[64];
  long val = -1, avail = -1;
  while (std::fscanf(f, "%63s %ld kB\n", key, &val) == 2)
    if (std::strcmp(key, "MemAvailable:") == 0) { avail = val; break; }
  std::fclose(f);
  return avail;
}

VolumeView::VolumeView(QWidget *parent) : QVTKOpenGLNativeWidget(parent) {
  setRenderWindow(renderWindow_);
  renderWindow_->AddRenderer(renderer_);
  renderer_->SetBackground(0.08, 0.08, 0.10);

  // grayscale transfer function: intensity -> gray, with opacity ramp
  color_->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
  color_->AddRGBPoint(255.0, 1.0, 1.0, 1.0);
  opacity_->AddPoint(0.0, 0.0);
  opacity_->AddPoint(64.0, 0.0);
  opacity_->AddPoint(255.0, 0.85);

  property_->SetColor(color_);
  property_->SetScalarOpacity(opacity_);
  property_->SetInterpolationTypeToLinear();
  property_->ShadeOff();

  mapper_->SetBlendModeToComposite();
  mapper_->SetRequestedRenderModeToDefault();  // GPU if supported, else CPU fallback
  mapper_->SetAutoAdjustSampleDistances(1);    // coarser sampling while interacting
  mapper_->SetUseJittering(1);                  // break up wood-grain ray-cast ringing
  renderer_->UseFXAAOn();                        // cheap edge antialiasing
  volume_->SetMapper(mapper_);
  volume_->SetProperty(property_);
  // volume actor is added to the renderer on first setVolume() (avoids rendering
  // a mapper with no input before any data is loaded)
}

void VolumeView::setVolume(const uint8_t *data, int dz, int dy, int dx) {
  image_->SetDimensions(dx, dy, dz);
  image_->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
  std::memcpy(image_->GetScalarPointer(), data, (size_t)dz * dy * dx);
  image_->Modified();

  mapper_->SetInputData(image_);
  if (!added_) { renderer_->AddVolume(volume_); added_ = true; }
  renderer_->ResetCamera();
  renderWindow_->Render();
}

QString VolumeView::gpuInfo() {
  makeCurrent();
  auto *f = QOpenGLContext::currentContext()->functions();
  auto str = [f](GLenum e) {
    const char *s = reinterpret_cast<const char *>(f->glGetString(e));
    return s ? QString(s) : QStringLiteral("?");
  };
  GLint max3d = 0, max2d = 0;
  f->glGetIntegerv(GL_MAX_3D_TEXTURE_SIZE, &max3d);
  f->glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max2d);

  // dedicated-VRAM query, if the driver exposes one
  long vram_total_kb = -1, vram_avail_kb = -1;
  QString memSource = QStringLiteral("integrated / shared system RAM");
  if (auto *ctx = context()) {
    if (ctx->hasExtension(QByteArrayLiteral("GL_NVX_gpu_memory_info"))) {
      GLint t = 0, a = 0;
      f->glGetIntegerv(GL_GPU_MEMORY_INFO_TOTAL_AVAILABLE_MEMORY_NVX, &t);
      f->glGetIntegerv(GL_GPU_MEMORY_INFO_CURRENT_AVAILABLE_VIDMEM_NVX, &a);
      vram_total_kb = t; vram_avail_kb = a;
      memSource = QStringLiteral("VRAM (NVX)");
    } else if (ctx->hasExtension(QByteArrayLiteral("GL_ATI_meminfo"))) {
      GLint info[4] = {0, 0, 0, 0};
      f->glGetIntegerv(GL_TEXTURE_FREE_MEMORY_ATI, info);
      vram_avail_kb = info[0];
      memSource = QStringLiteral("VRAM (ATI)");
    }
  }
  doneCurrent();

  // memory budget for a single u8 3D texture (60% headroom for fb/gradient/etc.)
  long budget_kb = (vram_avail_kb > 0) ? vram_avail_kb : sys_mem_available_kb();
  double budget_bytes = budget_kb > 0 ? (double)budget_kb * 1024.0 * 0.6 : 0.0;
  long mem_edge = budget_bytes > 0 ? (long)std::cbrt(budget_bytes) : 0;
  long max_edge = max3d > 0 ? std::min<long>(max3d, mem_edge ? mem_edge : max3d) : mem_edge;
  const char *limiter = (mem_edge && mem_edge < max3d) ? "memory" : "texture-dim cap";

  QString s;
  s += QStringLiteral("GPU:        %1  (%2)\n").arg(str(GL_RENDERER), str(GL_VENDOR));
  s += QStringLiteral("GL version: %1\n").arg(str(GL_VERSION));
  s += QStringLiteral("max 3D texture (per axis): %1     max 2D: %2\n").arg(max3d).arg(max2d);
  if (vram_total_kb > 0)
    s += QStringLiteral("VRAM: %1 MiB total, %2 MiB free (%3)\n")
             .arg(vram_total_kb / 1024).arg(vram_avail_kb / 1024).arg(memSource);
  else
    s += QStringLiteral("VRAM: %1 (avail ~%2 MiB)\n").arg(memSource).arg(budget_kb / 1024);
  s += QStringLiteral("→ max single u8 volume: ~%1³ voxels (limited by %2)\n")
           .arg(max_edge).arg(limiter);

  // VTK-specific: is GPU ray-cast actually supported on this context, or will
  // SmartVolumeMapper silently fall back to (slow) CPU ray casting?
  vtkNew<vtkGPUVolumeRayCastMapper> probe;
  bool gpuRC = probe->IsRenderSupported(renderWindow_, property_) != 0;
  s += QStringLiteral("GPU volume ray-cast supported: %1")
           .arg(gpuRC ? "YES" : "NO — would fall back to CPU ray cast");
  if (added_) {
    const char *mode = "other";
    switch (mapper_->GetLastUsedRenderMode()) {
      case vtkSmartVolumeMapper::GPURenderMode: mode = "GPU ray cast"; break;
      case vtkSmartVolumeMapper::RayCastRenderMode: mode = "CPU ray cast"; break;
      default: break;
    }
    s += QStringLiteral("\nlast render mode actually used: %1").arg(mode);
  }
  return s;
}
