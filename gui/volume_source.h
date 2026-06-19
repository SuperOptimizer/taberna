/* volume_source.h — the ONLY place the GUI touches matter-compressor.
 *
 * Opens an .mca archive (mmap + mc_open, the stable read-only reader side) and
 * serves dense u8 voxel regions / 2D slices to the viewers. Everything mc-specific
 * lives here so the ongoing matter-compressor rewrite only touches this file, and
 * so the Qt/VTK frontend stays backend-agnostic (swap this class for a different
 * VolumeSource and the viewers are unchanged).
 *
 * Layout: dense regions are z-major, x-fastest — out[(z*dy + y)*dx + x].
 */
#pragma once
#include <cstdint>
#include <string>
#include <vector>

class VolumeSource {
public:
  ~VolumeSource();
  bool open(const std::string &mca_path);   // false on failure
  void close();
  bool isOpen() const { return archive_ != nullptr; }

  int   nx() const { return nx_; }
  int   ny() const { return ny_; }
  int   nz() const { return nz_; }
  int   nlods() const { return nlods_; }
  float quality() const { return quality_; }

  // voxel dims at a LOD (each level halves; clamped >= 1)
  int nxAt(int lod) const { return dimAt(nx_, lod); }
  int nyAt(int lod) const { return dimAt(ny_, lod); }
  int nzAt(int lod) const { return dimAt(nz_, lod); }

  // Decode a dense axis-aligned region of `lod` into `out` (caller-sized
  // dz*dy*dx). Out-of-range / absent voxels are 0.
  void readRegion(int lod, int z0, int y0, int x0, int dz, int dy, int dx,
                  uint8_t *out) const;

  // 2D slice at `lod`. axis: 0=axial (const z), 1=coronal (const y),
  // 2=sagittal (const x). Fills `out` (w*h, row-major) and sets w,h.
  void readSlice(int axis, int index, int lod, std::vector<uint8_t> &out,
                 int &w, int &h) const;

private:
  static int dimAt(int d, int lod) { int v = d >> lod; return v < 1 ? 1 : v; }

  struct mc_archive *archive_ = nullptr;   // read-region handle (kept open)
  int            nx_ = 0, ny_ = 0, nz_ = 0, nlods_ = 0;
  float          quality_ = 0.0f;
};
