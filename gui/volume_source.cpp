/* volume_source.cpp — see volume_source.h. The only matter-compressor coupling. */
#include "volume_source.h"

#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

extern "C" {
#include "matter_compressor.h"
}

VolumeSource::~VolumeSource() { close(); }

bool VolumeSource::open(const std::string &path) {
  close();

  // 1) mmap the archive once and read header metadata via the read-only reader.
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0) return false;
  struct stat st;
  if (fstat(fd, &st) != 0 || st.st_size <= 0) { ::close(fd); return false; }
  size_t len = (size_t)st.st_size;
  void *map = mmap(nullptr, len, PROT_READ, MAP_PRIVATE, fd, 0);
  if (map == MAP_FAILED) { ::close(fd); return false; }

  mc_reader *r = mc_open((const uint8_t *)map, len);
  if (!r) { munmap(map, len); ::close(fd); return false; }
  mc_reader_dims(r, &nx_, &ny_, &nz_);
  nlods_ = mc_reader_nlods(r);
  quality_ = mc_reader_quality(r);
  mc_close(r);
  munmap(map, len);
  ::close(fd);

  // 2) open an archive handle for region decoding (dims/quality must match header).
  archive_ = mc_archive_open_dims(path.c_str(), nx_, ny_, nz_, quality_);
  return archive_ != nullptr;
}

void VolumeSource::close() {
  if (archive_) { mc_archive_close(archive_); archive_ = nullptr; }
  nx_ = ny_ = nz_ = nlods_ = 0;
  quality_ = 0.0f;
}

void VolumeSource::readRegion(int lod, int z0, int y0, int x0, int dz, int dy, int dx,
                              uint8_t *out) const {
  if (!archive_) { std::memset(out, 0, (size_t)dz * dy * dx); return; }
  // dense C-order z-major x-fastest: out[(z*dy + y)*dx + x]
  size_t sy = (size_t)dx;
  size_t sz = (size_t)dx * (size_t)dy;
  mc_archive_read_region(archive_, lod, z0, y0, x0, dz, dy, dx, out, sz, sy, 0);
}

void VolumeSource::readSlice(int axis, int index, int lod, std::vector<uint8_t> &out,
                             int &w, int &h) const {
  int z0 = 0, y0 = 0, x0 = 0;
  int dz = nzAt(lod), dy = nyAt(lod), dx = nxAt(lod);
  switch (axis) {
    case 0: z0 = index; dz = 1; w = dx; h = dy; break;  // axial   (const z)
    case 1: y0 = index; dy = 1; w = dx; h = dz; break;  // coronal (const y)
    default: x0 = index; dx = 1; w = dy; h = dz; break; // sagittal(const x)
  }
  out.assign((size_t)dz * dy * dx, 0);
  readRegion(lod, z0, y0, x0, dz, dy, dx, out.data());
  // the singleton axis collapses, leaving a dense w*h row-major image
}
