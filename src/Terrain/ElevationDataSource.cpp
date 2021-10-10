/*
 *  Copyright 2019-2021 Diligent Graphics LLC
 *  Copyright 2015-2019 Egor Yusov
 *  
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *  
 *      http://www.apache.org/licenses/LICENSE-2.0
 *  
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  In no event and under no legal theory, whether in tort (including negligence), 
 *  contract, or otherwise, unless required by applicable law (such as deliberate 
 *  and grossly negligent acts) or agreed to in writing, shall any Contributor be
 *  liable for any damages, including any direct, indirect, special, incidental, 
 *  or consequential damages of any character arising as a result of this License or 
 *  out of the use or inability to use the software (including but not limited to damages 
 *  for loss of goodwill, work stoppage, computer failure or malfunction, or any and 
 *  all other commercial damages or losses), even if such Contributor has been advised 
 *  of the possibility of such damages.
 */

// This file is derived from the open source project provided by Intel Corportaion that
// requires the following notice to be kept:
//--------------------------------------------------------------------------------------
// Copyright 2013 Intel Corporation
// All Rights Reserved
//
// Permission is granted to use, copy, distribute and prepare derivative works of this
// software for any purpose and without fee, provided, that the above copyright notice
// and this statement appear in all copies.  Intel makes no representations about the
// suitability of this software for any purpose.  THIS SOFTWARE IS PROVIDED "AS IS."
// INTEL SPECIFICALLY DISCLAIMS ALL WARRANTIES, EXPRESS OR IMPLIED, AND ALL LIABILITY,
// INCLUDING CONSEQUENTIAL AND OTHER INDIRECT DAMAGES, FOR THE USE OF THIS SOFTWARE,
// INCLUDING LIABILITY FOR INFRINGEMENT OF ANY PROPRIETARY RIGHTS, AND INCLUDING THE
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  Intel does not
// assume any responsibility for any errors which may appear in this software nor any
// responsibility to update it.
//--------------------------------------------------------------------------------------

#include <algorithm>
#include <cmath>

#include "BasicFileStream.hpp"
#include "DataBlobImpl.hpp"
#include "ElevationDataSource.hpp"
#include "FileWrapper.hpp"
#include "GraphicsAccessories.hpp"
#include "Image.h"
#include "TextureUtilities.h"

namespace Diligent
{

// Creates data source from the specified raw data file
ElevationDataSource::ElevationDataSource(const Char* src_file) : num_levels_(0), patch_size_(128), col_offset_(0), row_offset_(0)
{

#if 1 // read height map from file
  RefCntAutoPtr<Image> height_map_image;
  CreateImageFromFile(src_file, &height_map_image);

  const auto& img_info   = height_map_image->GetDesc();
  auto*       image_data = height_map_image->GetData();
  // Calculate minimal number of columns and rows in the form 2^n+1 that encompass the data to ensure the mask map repeateable
  dim_col_ = 1;
  dim_row_ = 1;
  while (dim_col_ + 1 < img_info.Width || dim_row_ + 1 < img_info.Height)
  {
    dim_col_ *= 2;
    dim_row_ *= 2;
  }

  num_levels_ = 1;
  while ((patch_size_ << (num_levels_ - 1)) < dim_col_ || (patch_size_ << (num_levels_ - 1)) < dim_row_) num_levels_++;

  dim_col_++;
  dim_row_++;
  stride_ = dim_col_;
  //stride_ = (num_cols_ + 1) & (-2); // -2 is 1111 1110 , to ensure num_cols_ is odd

  // Load the data
  height_map_.resize(dim_row_ * static_cast<size_t>(stride_));

  VERIFY(img_info.ComponentType == VT_UINT16 && img_info.NumComponents == 1, "Unexpected scanline size: 16-bit single-channel image is expected");
  auto* img_data = reinterpret_cast<Uint8*>(image_data->GetDataPtr());
  for (Uint32 row = 0; row < img_info.Height; row++, img_data += img_info.RowStride)
  { memcpy(&height_map_[row * static_cast<size_t>(stride_)], img_data, size_t{img_info.Width} * size_t{GetValueSize(img_info.ComponentType)}); }

  // Duplicate the last row and column, why?
  for (Uint32 iRow = 0; iRow < img_info.Height; iRow++)
    for (Uint32 iCol = img_info.Width; iCol < dim_col_; iCol++) GetHeight(iCol, iRow) = GetHeight((img_info.Width - 1), iRow);

  for (Uint32 iCol = 0; iCol < dim_col_; iCol++)
    for (Uint32 iRow = img_info.Height; iRow < dim_row_; iRow++) GetHeight(iCol, iRow) = GetHeight(iCol, img_info.Height - 1);

#else
  stride_   = 2048;
  num_rows_ = num_cols_ = 2048;
  num_levels_           = 5;
  height_map_.resize(num_rows_ * static_cast<size_t>(stride_));
  for (Uint32 j = 0; j < num_rows_; ++j)
  {
    for (Uint32 i = 0; i < num_cols_; ++i)
    {
      float x = (float)i / (float)num_rows_;
      float y = (float)j / (float)num_cols_;
      float A = 1;
      float F = 1;
      float h = 0;
      for (int Octave = 0; Octave < 8; ++Octave, A *= 0.7f, F *= 1.8f)
      {
        h += A * (sin(x * F) * sin(y * 1.5f * F) + 0.5f * (sin((x + y) * 1.3f * F) + cos((x * y) * 1.7f * F)));
        h = fabs(h - 0.5f);
      }
      h = fabs(h) * 32000.f;
      h = std::min(h, (float)std::numeric_limits<Uint16>::max());

      GetElevSample(i, j) = (Uint16)h;
    }
  }
#endif
  min_max_elevation_.Resize(num_levels_);

  // Calculate min/max elevations
  CalculateMinMaxElevations();
}

ElevationDataSource::~ElevationDataSource(void) {}

Uint16 ElevationDataSource::GetGlobalMinElevation() const { return min_max_elevation_[QuadTreeNodePosition()].first; }

Uint16 ElevationDataSource::GetGlobalMaxElevation() const { return min_max_elevation_[QuadTreeNodePosition()].second; }

int MirrorCoord(int coord, int dim)
{
  coord      = std::abs(coord);
  int period = coord / dim;
  coord      = coord % dim;
  if (period & 0x01) { coord = (dim - 1) - coord; }
  return coord;
}

inline Uint16& ElevationDataSource::GetHeight(size_t col, size_t row) { return height_map_[col + row * static_cast<size_t>(stride_)]; }

inline Uint16 ElevationDataSource::GetHeight(size_t col, size_t row) const { return height_map_[col + row * static_cast<size_t>(stride_)]; }

/*
¡¡¡¡step: distance of the data used for interpolated calculation
*/
float ElevationDataSource::GetInterpolateHeight(float col, float row, int step) const
{
  int   col_0 = static_cast<int>(floor(col / static_cast<float>(step))) * step;
  int   row_0 = static_cast<int>(floor(row / static_cast<float>(step))) * step;
  float w_col = (col - static_cast<float>(col_0)) / static_cast<float>(step);
  float w_row = (row - static_cast<float>(row_0)) / static_cast<float>(step);
  col_0 += col_offset_;
  row_0 += row_offset_;

  int col_1 = col_0 + step;
  int row_1 = row_0 + step;

  col_0 = MirrorCoord(col_0, dim_col_);
  col_1 = MirrorCoord(col_1, dim_col_);
  row_0 = MirrorCoord(row_0, dim_row_);
  row_1 = MirrorCoord(row_1, dim_row_);

  Uint16 h_00 = GetHeight(col_0, row_0);
  Uint16 h_10 = GetHeight(col_1, row_0);
  Uint16 h_01 = GetHeight(col_0, row_1);
  Uint16 h_11 = GetHeight(col_1, row_1);

  return (h_00 * w_col + h_10 * (1 - w_col)) * w_row + (h_01 * w_col + h_11 * (1 - w_col)) * (1 - w_row);
  /*return (h_00 * (1 - w_col) + h_10 * w_col) * (1 - w_row) + (h_01 * (1 - w_col) + h_11 * w_col) * w_row;*/
}

/*
pixel_size: terrain map unit, or scale for each pixel
*/
float3 ElevationDataSource::ComputeNormal(float col, float row, float pixel_size, float height_scale, int step) const
{
  float float_step = static_cast<float>(step);
  float height_1   = GetInterpolateHeight(col + float_step, row, step);
  float hieght_2   = GetInterpolateHeight(col - float_step, row, step);
  float height_3   = GetInterpolateHeight(col, row + float_step, step);
  float height_4   = GetInterpolateHeight(col, row - float_step, step);

  float3 grad;
  grad.x = hieght_2 - height_1;
  grad.y = height_4 - height_3;
  grad.z = float_step * pixel_size * 2.f;

  grad.x *= height_scale;
  grad.y *= height_scale;
  return normalize(grad);
}

void ElevationDataSource::RecomputePatchMinMaxElevations(const QuadTreeNodePosition& pos)
{
  if (pos.level == num_levels_ - 1)
  {
    std::pair<Uint16, Uint16>& now_patch_min_max_elev = min_max_elevation_[QuadTreeNodePosition(pos.horz_order, pos.vert_order, pos.level)];

    int iStartCol = pos.horz_order * patch_size_;
    int iStartRow = pos.vert_order * patch_size_;

    now_patch_min_max_elev.first = now_patch_min_max_elev.second = GetHeight(iStartCol, iStartRow);
    for (Uint32 iRow = iStartRow; iRow <= iStartRow + patch_size_; iRow++)
      for (Uint32 iCol = iStartCol; iCol <= iStartCol + patch_size_; iCol++)
      {
        Uint16 now_elev               = GetHeight(std::min(iCol, (Uint32)dim_col_ - 1), std::min(iRow, (Uint32)dim_row_ - 1));
        now_patch_min_max_elev.first  = std::min(now_patch_min_max_elev.first, now_elev);
        now_patch_min_max_elev.second = std::max(now_patch_min_max_elev.second, now_elev);
      }
  }
  else
  {
    std::pair<Uint16, Uint16>& now_patch_min_max_elev = min_max_elevation_[pos];

    std::pair<Uint16, Uint16>& lb_child_min_max_elev = min_max_elevation_[GetChildPosition(pos, 0)];
    std::pair<Uint16, Uint16>& rb_child_min_max_elev = min_max_elevation_[GetChildPosition(pos, 1)];
    std::pair<Uint16, Uint16>& lt_child_min_max_elev = min_max_elevation_[GetChildPosition(pos, 2)];
    std::pair<Uint16, Uint16>& rt_child_min_max_elev = min_max_elevation_[GetChildPosition(pos, 3)];

    now_patch_min_max_elev.first = std::min(lb_child_min_max_elev.first, rb_child_min_max_elev.first);
    now_patch_min_max_elev.first = std::min(now_patch_min_max_elev.first, lt_child_min_max_elev.first);
    now_patch_min_max_elev.first = std::min(now_patch_min_max_elev.first, rt_child_min_max_elev.first);

    now_patch_min_max_elev.second = std::max(lb_child_min_max_elev.second, rb_child_min_max_elev.second);
    now_patch_min_max_elev.second = std::max(now_patch_min_max_elev.second, lt_child_min_max_elev.second);
    now_patch_min_max_elev.second = std::max(now_patch_min_max_elev.second, rt_child_min_max_elev.second);
  }
}

// Calculates min/max elevations for the hierarchy
void ElevationDataSource::CalculateMinMaxElevations()
{
  // Calculate min/max elevations starting from the finest level
  for (HierarchyReverseIterator it(num_levels_); it.IsValid(); it.Next()) { RecomputePatchMinMaxElevations(it); }
}

void ElevationDataSource::GetDataPtr(const Uint16*& data_ptr, size_t& stride)
{
  data_ptr = &height_map_[0];
  stride   = stride_;
}

} // namespace Diligent
