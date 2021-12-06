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

#include "EarthHemisphere.hpp"

#include <algorithm>
#include <array>
#include <cfloat>
#include <type_traits>

namespace Diligent
{

#include "BasicStructures.fxh"
// clang-format off
// ToneMappingStructures.fxh must be included before EpipolarLightScatteringStructures.fxh
#include "ToneMappingStructures.fxh"
#include "EpipolarLightScatteringStructures.fxh"
// clang-format on

} // namespace Diligent

#include "CommonlyUsedStates.h"
#include "ElevationDataSource.hpp"
#include "GraphicsAccessories.hpp"
#include "GraphicsUtilities.h"
#include "MapHelper.hpp"
#include "ShaderMacroHelper.hpp"
#include "TextureUtilities.h"

namespace Diligent
{

struct HemisphereVertex
{
  float3 world_pos;
  float2 mask_uv; // mask coord for normal map?
  HemisphereVertex() : world_pos(0, 0, 0), mask_uv(0, 0) {}
};

template <typename E>
constexpr auto to_underlying(E e) noexcept
{
  return static_cast<std::underlying_type_t<E>>(e);
}

enum class QUAD_TRI_TYPE : Uint8
{
  UNDEFINED = 2,
  // 01      11 (row 1)
  //  *------*
  //  |   .' |
  //  | .'   |
  //  * -----*
  // 00     10  (row 0)
  DIAG_00_11 = 1, // where frist index is in row 1
  // 01      11 (row 1)
  //  *------*
  //  | '.   |
  //  |   '. |
  //  * -----*
  // 00      10 (row 0)
  DIAG_01_10 = 0, // where frist index is in row 0
};

template <typename IndexType, template <typename> class IndexGenerator>
class TriStrip
{
public:
  TriStrip(std::vector<IndexType>& indices, IndexGenerator<IndexType> index_generator) :
    quad_tri_type_(QUAD_TRI_TYPE::UNDEFINED), indices_(indices), index_generator_(index_generator)
  {}

  void AddStrip(IndexType base_index, IndexType init_col, IndexType init_row, IndexType num_cols, IndexType num_rows, QUAD_TRI_TYPE quad_tri_type)
  {
    VERIFY_EXPR(quad_tri_type == QUAD_TRI_TYPE::DIAG_00_11 || quad_tri_type == QUAD_TRI_TYPE::DIAG_01_10);
    IndexType first_index = base_index + index_generator_(init_col, init_row + static_cast<IndexType>(quad_tri_type));
    if (quad_tri_type_ != QUAD_TRI_TYPE::UNDEFINED)
    {
      // To move from one strip to another, we have to generate two degenerate triangles
      // by duplicating the last vertex in previous strip and the first vertex in new strip
      indices_.push_back(indices_.back());
      indices_.push_back(first_index);
    }

    if ((quad_tri_type_ != QUAD_TRI_TYPE::UNDEFINED && quad_tri_type_ != quad_tri_type) ||
        (quad_tri_type_ == QUAD_TRI_TYPE::UNDEFINED && quad_tri_type == QUAD_TRI_TYPE::DIAG_01_10))
    {
      // If triangulation orientation changes, or if start strip orientation is 01 to 10,
      // we also have to add one additional vertex to preserve winding order
      indices_.push_back(first_index);
    }
    quad_tri_type_ = quad_tri_type;

    for (IndexType row = 0; row < num_rows - 1; ++row)
    {
      for (IndexType col = 0; col < num_cols; ++col)
      {
        IndexType v00 = base_index + index_generator_(init_col + col, init_row + row);
        IndexType v01 = base_index + index_generator_(init_col + col, init_row + row + 1);
        if (quad_tri_type_ == QUAD_TRI_TYPE::DIAG_01_10)
        {
          if (col == 0 && row == 0) VERIFY_EXPR(first_index == v00);
          // 01      11
          //  *------*
          //  | '.   |
          //  |   '. |
          //  * -----*
          // 00      10
          indices_.push_back(v00);
          indices_.push_back(v01);
        }
        else if (quad_tri_type_ == QUAD_TRI_TYPE::DIAG_00_11)
        {
          if (col == 0 && row == 0) VERIFY_EXPR(first_index == v01);
          // 01      11
          //  *------*
          //  |   .' |
          //  | .'   |
          //  * -----*
          // 00      10
          indices_.push_back(v01);
          indices_.push_back(v00);
        }
        else
        {
          VERIFY_EXPR(false);
        }
      }

      if (row < num_rows - 2)
      {
        indices_.push_back(indices_.back());
        indices_.push_back(base_index + index_generator_(init_col, init_row + row + static_cast<IndexType>(quad_tri_type) + 1));
      }
    }
  }

private:
  QUAD_TRI_TYPE             quad_tri_type_;
  std::vector<IndexType>&   indices_;
  IndexGenerator<IndexType> index_generator_;
};

template <typename IndexType>
class StdIndexGenerator
{
public:
  StdIndexGenerator(IndexType dim) : dim_(dim) {}

  IndexType operator()(IndexType col, IndexType row, IndexType init = 0) { return col + row * dim_ + init; }

private:
  IndexType dim_;
};

typedef TriStrip<Uint32, StdIndexGenerator> StdTriStrip32;

void ComputeVertexHeight(HemisphereVertex& vertex, class ElevationDataSource* elev_data_src, float sampling_step, float height_scale)
{
  float3& world_pos = vertex.world_pos;

  float hm_col   = world_pos.x / sampling_step; // sampling_step is the size by which the height map represents
  float hm_row   = world_pos.z / sampling_step;
  float altitude = elev_data_src->GetInterpolateHeight(hm_col, hm_row);
  int   col_offset, row_offset;
  elev_data_src->GetOffsets(col_offset, row_offset);
  // +0.5f due to translate to orign£¿
  vertex.mask_uv.x = (hm_col + static_cast<float>(col_offset) + 0.5f) / static_cast<float>(elev_data_src->GetDimCol());
  vertex.mask_uv.y = (hm_row + static_cast<float>(row_offset) + 0.5f) / static_cast<float>(elev_data_src->GetDimRow());

  float3 sphere_normal = normalize(world_pos);
  world_pos += sphere_normal * altitude * height_scale;
}

template <typename IndexType>
class RingMeshBuilder
{
public:
  RingMeshBuilder(IRenderDevice* device, const std::vector<HemisphereVertex>& vertex_buffer, int grid_dim, std::vector<RingSectorMesh>& ring_meshes) :
    device_(device), ring_meshes_(ring_meshes), vertex_buffer_(vertex_buffer), grid_dim_(grid_dim)
  {}

  void CreateMesh(IndexType base_index, IndexType init_col, IndexType init_row, IndexType num_cols, IndexType num_rows, QUAD_TRI_TYPE quad_tri_type)
  {
    ring_meshes_.push_back(RingSectorMesh());
    RingSectorMesh& now_mesh = ring_meshes_.back();

    std::vector<IndexType> index_buffer;
    StdTriStrip32          tri_strip(index_buffer, StdIndexGenerator<IndexType>(grid_dim_));
    tri_strip.AddStrip(base_index, init_col, init_row, num_cols, num_rows, quad_tri_type);

    now_mesh.num_indices = index_buffer.size();

    // Prepare buffer description
    BufferDesc index_buffer_desc;
    index_buffer_desc.Name      = "Ring mesh index buffer";
    index_buffer_desc.Size      = index_buffer.size() * sizeof(index_buffer[0]);
    index_buffer_desc.BindFlags = BIND_INDEX_BUFFER;
    index_buffer_desc.Usage     = USAGE_IMMUTABLE;
    BufferData index_buffer_data;
    index_buffer_data.pData    = index_buffer.data();
    index_buffer_data.DataSize = index_buffer_desc.Size;
    // Create the buffer
    device_->CreateBuffer(index_buffer_desc, &index_buffer_data, &now_mesh.index_buffer);
    VERIFY(now_mesh.index_buffer, "Failed to create index buffer");

    // Compute bounding box
    auto& bound_box = now_mesh.bound_box;
    bound_box.Max   = float3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    bound_box.Min   = float3(+FLT_MAX, +FLT_MAX, +FLT_MAX);
    for (auto i = index_buffer.begin(); i != index_buffer.end(); ++i)
    {
      const auto& now_vert_pos = vertex_buffer_[*i].world_pos;

      bound_box.Min = std::min(bound_box.Min, now_vert_pos);
      bound_box.Max = std::max(bound_box.Max, now_vert_pos);
    }
  }

private:
  RefCntAutoPtr<IRenderDevice>         device_;
  std::vector<RingSectorMesh>&         ring_meshes_;
  const std::vector<HemisphereVertex>& vertex_buffer_; // vertex buffer
  const IndexType                      grid_dim_;
};

template <typename IndexType>
void GenerateSphereGeometry(IRenderDevice*                 device,
                            const float                    earth_radius,
                            IndexType                      grid_dim, // grid dimension
                            const size_t                   num_rings,
                            class ElevationDataSource*     elev_data_src,
                            float                          sampling_step,
                            float                          height_scale,
                            std::vector<HemisphereVertex>& vert_buff,
                            std::vector<RingSectorMesh>&   sphere_meshes)
{
  // Grid dimension of projection sphere, due to symmetry, it should be 2K+1, but why 4K+1
  if ((grid_dim - 1) % 4 != 0)
  {
    grid_dim = RenderingParams().ring_dim;
    UNEXPECTED("Grid dimension must be 4k+1. Defaulting to ", grid_dim);
  }
  const IndexType           grid_midst = (grid_dim - 1) / 2;
  const IndexType           grid_quart = (grid_dim - 1) / 4;
  StdIndexGenerator<Uint32> std_index_generator(grid_dim);
  RingMeshBuilder<Uint32>   ring_mesh_builder(device, vert_buff, grid_dim, sphere_meshes);

  size_t start_ring = 0;
  vert_buff.reserve(num_rings * grid_dim * grid_dim);
  for (size_t ring = start_ring; ring < num_rings; ++ring)
  {
    size_t now_index_init = vert_buff.size();
    // Why not resize once directly? the actual used size might < grid_dim * grid_dim ?
    vert_buff.resize(vert_buff.size() + static_cast<size_t>(grid_dim) * static_cast<size_t>(grid_dim));
    // What is grid scale? grid_scale = 2^-(num_rings-1-i_ring)
    float grid_scale = 1.f / static_cast<float>(1 << (num_rings - 1 - ring));
    // Fill vertex buffer
    for (size_t row = 0; row < grid_dim; ++row)
      for (size_t col = 0; col < grid_dim; ++col)
      {
        HemisphereVertex& now_vert = vert_buff[std_index_generator(col, row, now_index_init)];
        float3&           pos      = now_vert.world_pos;
        // Relative postion in a grid_dim * grid_dim grid
        pos.x = static_cast<float>(col) / static_cast<float>(grid_dim - 1);
        pos.z = static_cast<float>(row) / static_cast<float>(grid_dim - 1);
        pos.x = pos.x * 2 - 1; // Translate to origin, i.e. the sphere center
        pos.z = pos.z * 2 - 1;
        pos.y = 0;

        // Direction scale: projection from a square to a circle along
        // min/max direction
        float direction_scale = 1;
        if (pos.x != 0 || pos.z != 0)
        {
          float magnitude_x = fabs(pos.x);
          float magnitude_y = fabs(pos.z);
          float max_leg     = std::max(magnitude_x, magnitude_y);
          float min_leg     = std::min(magnitude_x, magnitude_y);
          float tan         = min_leg / max_leg;
          direction_scale   = 1 / sqrt(1 + tan * tan); // Cos value along min/max direction
        }

        pos.x *= direction_scale * grid_scale;
        pos.z *= direction_scale * grid_scale;
        pos.y = sqrt(std::max(0.f, 1.f - (pos.x * pos.x + pos.z * pos.z))); // Assure in an unit sphere

        pos.x *= earth_radius;
        pos.z *= earth_radius;
        pos.y *= earth_radius;

        ComputeVertexHeight(now_vert, elev_data_src, sampling_step, height_scale);
        pos.y -= earth_radius; // Translate the top to the origin
      }

    // Align vertices on the outer boundary, it can stil run if we skip this section
    if (ring < num_rings - 1)
    {
      for (size_t i = 1; i < grid_dim - 1; i += 2)
      {
        // Top & bottom boundaries
        for (size_t row = 0; row < grid_dim; row += grid_dim - 1)
        {
          const auto& V0 = vert_buff[std_index_generator(i - 1, row, now_index_init)].world_pos;
          auto&       V1 = vert_buff[std_index_generator(i + 0, row, now_index_init)].world_pos;
          const auto& V2 = vert_buff[std_index_generator(i + 1, row, now_index_init)].world_pos;
          V1             = (V0 + V2) / 2.f;
        }

        // Left & right boundaries
        for (size_t col = 0; col < grid_dim; col += grid_dim - 1)
        {
          const auto& V0 = vert_buff[std_index_generator(col, i - 1, now_index_init)].world_pos;
          auto&       V1 = vert_buff[std_index_generator(col, i + 0, now_index_init)].world_pos;
          const auto& V2 = vert_buff[std_index_generator(col, i + 1, now_index_init)].world_pos;
          V1             = (V0 + V2) / 2.f;
        }
      }
    }

    // Generate indices for the current ring
    if (ring == 0)
    {
      //ring_mesh_builder.CreateMesh(now_index_init, 0         , 0         , grid_dim +1, grid_dim + 1, QUAD_TRI_TYPE::DIAG_00_11);?
      // clang-format off
			ring_mesh_builder.CreateMesh(now_index_init, 0         , 0         , grid_midst + 1, grid_midst + 1, QUAD_TRI_TYPE::DIAG_00_11);
			ring_mesh_builder.CreateMesh(now_index_init, grid_midst, 0         , grid_midst + 1, grid_midst + 1, QUAD_TRI_TYPE::DIAG_01_10);
			ring_mesh_builder.CreateMesh(now_index_init, 0         , grid_midst, grid_midst + 1, grid_midst + 1, QUAD_TRI_TYPE::DIAG_01_10);
			ring_mesh_builder.CreateMesh(now_index_init, grid_midst, grid_midst, grid_midst + 1, grid_midst + 1, QUAD_TRI_TYPE::DIAG_00_11);
      // clang-format on
    }
    else
    {
      // clang-format off
			ring_mesh_builder.CreateMesh(now_index_init, 0         , 0, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);
			ring_mesh_builder.CreateMesh(now_index_init, grid_quart, 0, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);
      //ring_mesh_builder.CreateMesh(now_index_init, 0, 0, grid_midst + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);

			ring_mesh_builder.CreateMesh(now_index_init, grid_midst    , 0, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);
			ring_mesh_builder.CreateMesh(now_index_init, grid_quart * 3, 0, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);

			ring_mesh_builder.CreateMesh(now_index_init, 0, grid_quart, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);
			ring_mesh_builder.CreateMesh(now_index_init, 0, grid_midst, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);

			ring_mesh_builder.CreateMesh(now_index_init, grid_quart * 3, grid_quart, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);
			ring_mesh_builder.CreateMesh(now_index_init, grid_quart * 3, grid_midst, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);

			ring_mesh_builder.CreateMesh(now_index_init, 0         , grid_quart * 3, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);
			ring_mesh_builder.CreateMesh(now_index_init, grid_quart, grid_quart * 3, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);

			ring_mesh_builder.CreateMesh(now_index_init, grid_midst    , grid_quart * 3, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);
			ring_mesh_builder.CreateMesh(now_index_init, grid_quart * 3, grid_quart * 3, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);
      // clang-format on
    }
  }

  // We do not need per-vertex normals as we use normal map to shade terrain
  // Sphere tangent vertex are computed in the shader
#if 0
		// Compute normals
	const float3* pV0 = nullptr;
	const float3* pV1 = &vertex_buffer[index_buffer[0]].f3WorldPos;
	const float3* pV2 = &vertex_buffer[index_buffer[1]].f3WorldPos;
	float fSign = +1;
	for (Uint32 i = 2; i < m_uiIndicesInIndBuff; ++i)
	{
		fSign = -fSign;
		pV0 = pV1;
		pV1 = pV2;
		pV2 = &vertex_buffer[index_buffer[i]].f3WorldPos;
		float3 Rib0 = *pV0 - *pV1;
		float3 Rib1 = *pV1 - *pV2;
		float3 TriN;
		D3DXVec3Cross(&TriN, &Rib0, &Rib1);
		float fLength = D3DXVec3Length(&TriN);
		if (fLength > 0.1)
		{
			TriN /= fLength * fSign;
			for (int i = -2; i <= 0; ++i)
				vertex_buffer[index_buffer[i + i]].f3Normal += TriN;
		}
	}
	for (auto VBIt = vertex_buffer.begin(); VBIt != vertex_buffer.end(); ++VBIt)
	{
		float fLength = D3DXVec3Length(&VBIt->f3Normal);
		if (fLength > 1)
			VBIt->f3Normal /= fLength;
	}

	// Adjust normals on boundaries
	for (int iRing = iStartRing; iRing < iNumRings - 1; ++iRing)
	{
		int iCurrGridStart = (iRing - iStartRing) * iGridDimension * iGridDimension;
		int iNextGridStart = (iRing - iStartRing + 1) * iGridDimension * iGridDimension;
		for (int i = 0; i < iGridDimension; i += 2)
		{
			for (int Bnd = 0; Bnd < 2; ++Bnd)
			{
				const int CurrGridOffsets[] = { 0, iGridDimension - 1 };
				const int NextGridPffsets[] = { iGridQuart, iGridQuart * 3 };
				// Left and right boundaries
				{
					auto& CurrGridN = vertex_buffer[iCurrGridStart + CurrGridOffsets[Bnd] + i * iGridDimension].f3Normal;
					auto& NextGridN = vertex_buffer[iNextGridStart + NextGridPffsets[Bnd] + (iGridQuart + i / 2) * iGridDimension].f3Normal;
					auto NewN = CurrGridN + NextGridN;
					D3DXVec3Normalize(&NewN, &NewN);
					CurrGridN = NextGridN = NewN;
					if (i > 1)
					{
						auto& PrevCurrGridN = vertex_buffer[iCurrGridStart + CurrGridOffsets[Bnd] + (i - 2) * iGridDimension].f3Normal;
						auto MiddleN = PrevCurrGridN + NewN;
						D3DXVec3Normalize(&vertex_buffer[iCurrGridStart + CurrGridOffsets[Bnd] + (i - 1) * iGridDimension].f3Normal, &MiddleN);
					}
				}

				// Bottom and top boundaries
				{
					auto& CurrGridN = vertex_buffer[iCurrGridStart + i + CurrGridOffsets[Bnd] * iGridDimension].f3Normal;
					auto& NextGridN = vertex_buffer[iNextGridStart + (iGridQuart + i / 2) + NextGridPffsets[Bnd] * iGridDimension].f3Normal;
					auto NewN = CurrGridN + NextGridN;
					D3DXVec3Normalize(&NewN, &NewN);
					CurrGridN = NextGridN = NewN;
					if (i > 1)
					{
						auto& PrevCurrGridN = vertex_buffer[iCurrGridStart + (i - 2) + CurrGridOffsets[Bnd] * iGridDimension].f3Normal;
						auto MiddleN = PrevCurrGridN + NewN;
						D3DXVec3Normalize(&vertex_buffer[iCurrGridStart + (i - 1) + CurrGridOffsets[Bnd] * iGridDimension].f3Normal, &MiddleN);
					}
				}
			}
		}
	}
#endif
}

void EarthHemsiphere::RenderNormalMap(IRenderDevice*  device,
                                      IDeviceContext* context,
                                      const Uint16*   height_map_data,
                                      size_t          height_map_pitch,
                                      size_t          height_map_dim,
                                      ITexture*       normal_map_texture)
{
  TextureDesc height_map_desc;
  height_map_desc.Name      = "Height map texture";
  height_map_desc.Type      = RESOURCE_DIM_TEX_2D;
  height_map_desc.Width     = height_map_dim;
  height_map_desc.Height    = height_map_dim;
  height_map_desc.Format    = TEX_FORMAT_R16_UINT;
  height_map_desc.Usage     = USAGE_IMMUTABLE;
  height_map_desc.BindFlags = BIND_SHADER_RESOURCE;
  height_map_desc.MipLevels = ComputeMipLevelsCount(height_map_desc.Width, height_map_desc.Height);

  std::vector<Uint16> coarse_mip_levels;
  coarse_mip_levels.resize(height_map_dim / 2 * height_map_dim);

  std::vector<TextureSubResData> mip_map_data(height_map_desc.MipLevels);
  mip_map_data[0].pData        = height_map_data;
  mip_map_data[0].Stride       = (Uint32)height_map_pitch * sizeof(height_map_data[0]);
  const Uint16* finer_mip_data = height_map_data;
  Uint16*       now_mip_data   = &coarse_mip_levels[0];
  size_t        FinerMipPitch  = height_map_pitch;
  size_t        now_mip_pitch  = height_map_dim / 2; //All mip levels are stacked on top of each other with the same stride
  for (size_t mip_level = 1; mip_level < height_map_desc.MipLevels; ++mip_level)
  {
    auto mip_width  = height_map_desc.Width >> mip_level;
    auto mip_height = height_map_desc.Height >> mip_level;
    for (size_t row = 0; row < mip_height; ++row)
    {
      for (size_t col = 0; col < mip_width; ++col)
      {
        int average_height = 0;
        for (size_t i = 0; i < 2; ++i)
          for (size_t j = 0; j < 2; ++j) average_height += finer_mip_data[(col * 2 + i) + (row * 2 + j) * size_t{FinerMipPitch}];
        now_mip_data[col + row * now_mip_pitch] = (Uint16)(average_height >> 2);
      }
    }

    mip_map_data[mip_level].pData  = now_mip_data;
    mip_map_data[mip_level].Stride = (Uint32)now_mip_pitch * sizeof(*now_mip_data);
    finer_mip_data                 = now_mip_data;
    FinerMipPitch                  = now_mip_pitch;
    now_mip_data += mip_height * now_mip_pitch;
    now_mip_pitch = height_map_dim / 2;
  }

  RefCntAutoPtr<ITexture> height_map_texture;
  TextureData             height_map_texture_data;
  height_map_texture_data.pSubResources   = mip_map_data.data();
  height_map_texture_data.NumSubresources = (Uint32)mip_map_data.size();
  device->CreateTexture(height_map_desc, &height_map_texture_data, &height_map_texture);
  VERIFY(height_map_texture, "Failed to create height map texture");

  resource_mapping_->AddResource("g_tex2DElevationMap", height_map_texture->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE), true);

  RefCntAutoPtr<IBuffer> normal_generation_attribs_cbuffer;
  CreateUniformBuffer(device, sizeof(NMGenerationAttribs), "NM Generation Attribs CB", &normal_generation_attribs_cbuffer); //NM(normal)

  resource_mapping_->AddResource("cbNMGenerationAttribs", normal_generation_attribs_cbuffer, true);

  RefCntAutoPtr<IShaderSourceInputStreamFactory> shader_source_factory;
  device_->GetEngineFactory()->CreateDefaultShaderSourceStreamFactory("shaders\\;shaders\\terrain", &shader_source_factory);

  ShaderCreateInfo shader_cretae_info;
  shader_cretae_info.pShaderSourceStreamFactory = shader_source_factory;
  shader_cretae_info.FilePath                   = "ScreenSizeQuadVS.fx";
  shader_cretae_info.EntryPoint                 = "GenerateScreenSizeQuadVS";
  shader_cretae_info.SourceLanguage             = SHADER_SOURCE_LANGUAGE_HLSL;
  shader_cretae_info.UseCombinedTextureSamplers = true;
  shader_cretae_info.Desc.ShaderType            = SHADER_TYPE_VERTEX;
  shader_cretae_info.Desc.Name                  = "GenerateScreenSizeQuadVS";
  RefCntAutoPtr<IShader> screen_size_quad_vs;
  device->CreateShader(shader_cretae_info, &screen_size_quad_vs);

  shader_cretae_info.FilePath        = "GenerateNormalMapPS.fx";
  shader_cretae_info.EntryPoint      = "GenerateNormalMapPS";
  shader_cretae_info.Desc.ShaderType = SHADER_TYPE_PIXEL;
  shader_cretae_info.Desc.Name       = "GenerateNormalMapPS";

  RefCntAutoPtr<IShader> generate_normal_map_ps;
  device->CreateShader(shader_cretae_info, &generate_normal_map_ps);

  GraphicsPipelineStateCreateInfo pso_create_info;
  pso_create_info.pVS = screen_size_quad_vs;
  pso_create_info.pPS = generate_normal_map_ps;

  PipelineStateDesc& pso_desc                 = pso_create_info.PSODesc;
  pso_desc.ResourceLayout.DefaultVariableType = SHADER_RESOURCE_VARIABLE_TYPE_STATIC;
  // clang-format off
  ShaderResourceVariableDesc shader_vars[] = 
  {
    {SHADER_TYPE_PIXEL, "g_tex2DElevationMap",   SHADER_RESOURCE_VARIABLE_TYPE_STATIC},
    {SHADER_TYPE_PIXEL, "cbNMGenerationAttribs", SHADER_RESOURCE_VARIABLE_TYPE_STATIC}
  };
  // clang-format on
  pso_desc.ResourceLayout.NumVariables = _countof(shader_vars);
  pso_desc.ResourceLayout.Variables    = shader_vars;
  pso_desc.Name                        = "Render Normal Map";

  GraphicsPipelineDesc& graphics_pipline_desc                = pso_create_info.GraphicsPipeline;
  graphics_pipline_desc.DepthStencilDesc.DepthEnable         = false;
  graphics_pipline_desc.DepthStencilDesc.DepthWriteEnable    = false;
  graphics_pipline_desc.RasterizerDesc.FillMode              = FILL_MODE_SOLID;
  graphics_pipline_desc.RasterizerDesc.CullMode              = CULL_MODE_NONE;
  graphics_pipline_desc.RasterizerDesc.FrontCounterClockwise = true;
  graphics_pipline_desc.NumRenderTargets                     = 1;
  graphics_pipline_desc.RTVFormats[0]                        = TEX_FORMAT_RG8_UNORM;
  graphics_pipline_desc.PrimitiveTopology                    = PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP;
  RefCntAutoPtr<IPipelineState> render_normal_map_pso;
  device->CreateGraphicsPipelineState(pso_create_info, &render_normal_map_pso);
  render_normal_map_pso->BindStaticResources(SHADER_TYPE_VERTEX | SHADER_TYPE_PIXEL, resource_mapping_, BIND_SHADER_RESOURCES_VERIFY_ALL_RESOLVED);

  RefCntAutoPtr<IShaderResourceBinding> render_normal_map_srb;
  render_normal_map_pso->CreateShaderResourceBinding(&render_normal_map_srb, true);

  context->SetPipelineState(render_normal_map_pso);
  context->CommitShaderResources(render_normal_map_srb, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);

  const auto& normal_map_desc = normal_map_texture->GetDesc();
  for (Uint32 mip_level = 0; mip_level < normal_map_desc.MipLevels; ++mip_level)
  {
    TextureViewDesc texture_view_desc;
    texture_view_desc.ViewType        = TEXTURE_VIEW_RENDER_TARGET;
    texture_view_desc.MostDetailedMip = mip_level;
    RefCntAutoPtr<ITextureView> normal_map_rtv;
    normal_map_texture->CreateView(texture_view_desc, &normal_map_rtv);

    ITextureView* rtv_array[] = {normal_map_rtv};
    context->SetRenderTargets(_countof(rtv_array), rtv_array, nullptr, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);

    {
      MapHelper<NMGenerationAttribs> normal_generation_attribs(context, normal_generation_attribs_cbuffer, MAP_WRITE, MAP_FLAG_DISCARD);
      normal_generation_attribs->height_scale            = params_.terrain_attribs.height_scale;
      normal_generation_attribs->sample_spacing_interval = params_.terrain_attribs.elevation_sampling_interval;
      normal_generation_attribs->mip_level               = static_cast<int>(mip_level);
    }

    DrawAttribs draw_attrs(4, DRAW_FLAG_VERIFY_ALL);
    context->Draw(draw_attrs);
  }

  // Remove elevation map from resource mapping to release the resource
  resource_mapping_->RemoveResourceByName("g_tex2DElevationMap");
}

void EarthHemsiphere::Create(class ElevationDataSource* data_source,
                             const RenderingParams&     params,
                             IRenderDevice*             device,
                             IDeviceContext*            context,
                             const Char*                material_mask_path,
                             const Char*                tile_diffuse_path[],
                             const Char*                tile_normal_path[],
                             IBuffer*                   camera_attribs_buffer,
                             IBuffer*                   light_attribs_buffer,
                             IBuffer*                   media_scattering_params)
{
  params_ = params;
  device_ = device;

  const Uint16* height_map_data;
  size_t        height_map_pitch;
  data_source->GetDataPtr(height_map_data, height_map_pitch);
  Uint32 height_map_dim = data_source->GetDimCol();
  VERIFY_EXPR(height_map_dim == data_source->GetDimRow());

  TextureDesc normal_map_texture_desc;
  normal_map_texture_desc.Name      = "Normal map texture";
  normal_map_texture_desc.Type      = RESOURCE_DIM_TEX_2D;
  normal_map_texture_desc.Width     = height_map_dim;
  normal_map_texture_desc.Height    = height_map_dim;
  normal_map_texture_desc.Format    = TEX_FORMAT_RG8_UNORM;
  normal_map_texture_desc.Usage     = USAGE_DEFAULT;
  normal_map_texture_desc.BindFlags = BIND_SHADER_RESOURCE | BIND_RENDER_TARGET;
  normal_map_texture_desc.MipLevels = 0;

  RefCntAutoPtr<ITexture> normal_map_texture;
  device->CreateTexture(normal_map_texture_desc, nullptr, &normal_map_texture);
  normal_map_srv = normal_map_texture->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE);

  CreateUniformBuffer(device, sizeof(TerrainAttribs), "Terrain Attribs CB", &terrain_attribs_buffer_);

  ResourceMappingDesc res_mapping_desc;
  // clang-format off
  ResourceMappingEntry res_mapping_entries[] = 
  { 
    { "cbCameraAttribs", camera_attribs_buffer}, 
    { "cbTerrainAttribs", terrain_attribs_buffer_}, 
    { "cbLightAttribs", light_attribs_buffer}, 
    { "g_tex2DNormalMap", normal_map_texture->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE) }, 
    { "cbParticipatingMediaScatteringParams", media_scattering_params },
    {} 
  };
  // clang-format on
  res_mapping_desc.pEntries = res_mapping_entries;
  device->CreateResourceMapping(res_mapping_desc, &resource_mapping_);

  RefCntAutoPtr<ITexture> material_mask_texture;
  CreateTextureFromFile(material_mask_path, TextureLoadInfo(), device, &material_mask_texture);
  auto material_mask_srv = material_mask_texture->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE);
  resource_mapping_->AddResource("g_tex2DMtrlMap", material_mask_srv, true);

  // Load tiles
  IDeviceObject*          tile_diffuse_srv[NUM_TILE_TEXTURES] = {};
  RefCntAutoPtr<ITexture> tile_diffuse_texture[NUM_TILE_TEXTURES];
  IDeviceObject*          tile_normal_srv[NUM_TILE_TEXTURES] = {};
  RefCntAutoPtr<ITexture> tile_normal_texture[NUM_TILE_TEXTURES];
  for (size_t i = 0; i < static_cast<size_t>(NUM_TILE_TEXTURES); i++)
  {
    {
      TextureLoadInfo diff_tex_load_info;
      diff_tex_load_info.IsSRGB = false;
      CreateTextureFromFile(tile_diffuse_path[i], diff_tex_load_info, device, &tile_diffuse_texture[i]);
      tile_diffuse_srv[i] = tile_diffuse_texture[i]->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE);
    }

    {
      CreateTextureFromFile(tile_normal_path[i], TextureLoadInfo(), device, &tile_normal_texture[i]);
      tile_normal_srv[i] = tile_normal_texture[i]->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE);
    }
  }
  resource_mapping_->AddResourceArray("g_tex2DTileDiffuse", 0, tile_diffuse_srv, NUM_TILE_TEXTURES, true);
  resource_mapping_->AddResourceArray("g_tex2DTileNM", 0, tile_normal_srv, NUM_TILE_TEXTURES, true);

  device_->CreateSampler(Sam_ComparsionLinearClamp, &comparison_sampler_);

  RenderNormalMap(device, context, height_map_data, height_map_pitch, height_map_dim, normal_map_texture);

  RefCntAutoPtr<IShaderSourceInputStreamFactory> shader_source_factory;
  device_->GetEngineFactory()->CreateDefaultShaderSourceStreamFactory("shaders;shaders\\terrain;", &shader_source_factory);

  {
    ShaderCreateInfo shader_create_info;
    shader_create_info.pShaderSourceStreamFactory = shader_source_factory;
    shader_create_info.FilePath                   = "HemisphereVS.fx";
    shader_create_info.EntryPoint                 = "HemisphereVS";
    shader_create_info.SourceLanguage             = SHADER_SOURCE_LANGUAGE_HLSL;
    shader_create_info.UseCombinedTextureSamplers = true;
    shader_create_info.Desc.ShaderType            = SHADER_TYPE_VERTEX;
    shader_create_info.Desc.Name                  = "HemisphereVS";

    device->CreateShader(shader_create_info, &hemisphere_vs_);
  }

  {
    ShaderCreateInfo shader_create_info;
    shader_create_info.pShaderSourceStreamFactory = shader_source_factory;
    shader_create_info.FilePath                   = "HemisphereZOnlyVS.fx";
    shader_create_info.EntryPoint                 = "HemisphereZOnlyVS";
    shader_create_info.SourceLanguage             = SHADER_SOURCE_LANGUAGE_HLSL;
    shader_create_info.UseCombinedTextureSamplers = true;
    shader_create_info.Desc.ShaderType            = SHADER_TYPE_VERTEX;
    shader_create_info.Desc.Name                  = "HemisphereZOnlyVS";
    RefCntAutoPtr<IShader> hemisphere_z_only_vs_;
    device->CreateShader(shader_create_info, &hemisphere_z_only_vs_);

    GraphicsPipelineStateCreateInfo pso_create_info;
    PipelineStateDesc&              pso_desc = pso_create_info.PSODesc;

    pso_desc.Name                                         = "Render Hemisphere Z Only";
    pso_desc.ResourceLayout.DefaultVariableType           = SHADER_RESOURCE_VARIABLE_TYPE_STATIC;
    auto& GraphicsPipeline                                = pso_create_info.GraphicsPipeline;
    GraphicsPipeline.DepthStencilDesc                     = DSS_Default;
    GraphicsPipeline.RasterizerDesc.FillMode              = FILL_MODE_SOLID;
    GraphicsPipeline.RasterizerDesc.CullMode              = CULL_MODE_BACK;
    GraphicsPipeline.RasterizerDesc.DepthClipEnable       = False;
    GraphicsPipeline.RasterizerDesc.FrontCounterClockwise = True;
    // clang-format off
    LayoutElement layout_inputs[] =
    {
      {0, 0, 3, VT_FLOAT32, False, 0, (3+2)*4}
    };
    // clang-format on
    GraphicsPipeline.InputLayout.LayoutElements = layout_inputs;
    GraphicsPipeline.InputLayout.NumElements    = _countof(layout_inputs);
    GraphicsPipeline.DSVFormat                  = params_.ShadowMapFormat;
    GraphicsPipeline.PrimitiveTopology          = PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP;
    pso_create_info.pVS                         = hemisphere_z_only_vs_;
    device->CreateGraphicsPipelineState(pso_create_info, &hemisphere_z_only_pso_);
    hemisphere_z_only_pso_->BindStaticResources(SHADER_TYPE_VERTEX | SHADER_TYPE_PIXEL, resource_mapping_, BIND_SHADER_RESOURCES_VERIFY_ALL_RESOLVED);
    hemisphere_z_only_pso_->CreateShaderResourceBinding(&hemisphere_z_only_srb_, true);
  }

  std::vector<HemisphereVertex> vertex_buffer;
  GenerateSphereGeometry<Uint32>(device, Diligent::AirScatteringAttribs().fEarthRadius, params_.ring_dim, params_.num_rings, data_source,
                                 params_.terrain_attribs.elevation_sampling_interval, params_.terrain_attribs.height_scale, vertex_buffer,
                                 sphere_meshes_);

  BufferDesc vertex_buffer_desc;
  vertex_buffer_desc.Name      = "Hemisphere vertex buffer";
  vertex_buffer_desc.Size      = static_cast<Uint64>(vertex_buffer.size() * sizeof(vertex_buffer[0]));
  vertex_buffer_desc.Usage     = USAGE_IMMUTABLE;
  vertex_buffer_desc.BindFlags = BIND_VERTEX_BUFFER;
  BufferData vertex_buffer_init_data;
  vertex_buffer_init_data.pData    = vertex_buffer.data();
  vertex_buffer_init_data.DataSize = vertex_buffer_desc.Size;
  device->CreateBuffer(vertex_buffer_desc, &vertex_buffer_init_data, &vertex_buffer_);
  VERIFY(vertex_buffer_, "Failed to create vertex_buffer");
}

void EarthHemsiphere::Render(IDeviceContext*        context,
                             const RenderingParams& new_params,
                             const float4x4&        camera_view_proj_matrix,
                             ITextureView*          shadow_map_srv,
                             ITextureView*          pre_computeD_net_density_srv,
                             ITextureView*          ambient_skylight_srv,
                             bool                   z_only_pass)
{
  // clang-format off
  if (params_.num_shadow_cascades != new_params.num_shadow_cascades || 
      params_.best_cascade_search != new_params.best_cascade_search ||
      params_.filter_across_shadow_cascades != new_params.filter_across_shadow_cascades ||
      params_.fixed_shadow_filter_size != new_params.fixed_shadow_filter_size ||
      params_.DstRTVFormat != new_params.DstRTVFormat)
  // clang-format on
  {
    hemisphere_pso_.Release();
    hemisphere_srb_.Release();
  }

  params_ = new_params;

#if 0
    if( GetAsyncKeyState(VK_F9) )
    {
        m_RenderEarthHemisphereTech.Release();
    }
#endif

  if (!hemisphere_pso_)
  {
    ShaderCreateInfo shader_create_info;
    shader_create_info.FilePath                   = "HemispherePS.fx";
    shader_create_info.EntryPoint                 = "HemispherePS";
    shader_create_info.Desc.ShaderType            = SHADER_TYPE_PIXEL;
    shader_create_info.Desc.Name                  = "HemispherePS";
    shader_create_info.UseCombinedTextureSamplers = true;
    shader_create_info.SourceLanguage             = SHADER_SOURCE_LANGUAGE_HLSL;
    RefCntAutoPtr<IShaderSourceInputStreamFactory> shader_source_factory;
    device_->GetEngineFactory()->CreateDefaultShaderSourceStreamFactory("shaders;shaders\\terrain;", &shader_source_factory);
    shader_create_info.pShaderSourceStreamFactory = shader_source_factory;

    std::array<ImmutableSamplerDesc, 5> imtbl_samplers = {};

    imtbl_samplers[0].ShaderStages         = SHADER_TYPE_PIXEL;
    imtbl_samplers[0].SamplerOrTextureName = "g_tex2DTileDiffuse";
    imtbl_samplers[0].Desc.AddressU        = TEXTURE_ADDRESS_WRAP;
    imtbl_samplers[0].Desc.AddressV        = TEXTURE_ADDRESS_WRAP;
    imtbl_samplers[0].Desc.AddressW        = TEXTURE_ADDRESS_WRAP;

    imtbl_samplers[1].ShaderStages         = SHADER_TYPE_PIXEL;
    imtbl_samplers[1].SamplerOrTextureName = "g_tex2DTileNM";
    imtbl_samplers[1].Desc                 = imtbl_samplers[0].Desc;

    imtbl_samplers[2].ShaderStages         = SHADER_TYPE_PIXEL;
    imtbl_samplers[2].SamplerOrTextureName = "g_tex2DNormalMap";
    imtbl_samplers[2].Desc.AddressU        = TEXTURE_ADDRESS_MIRROR;
    imtbl_samplers[2].Desc.AddressV        = TEXTURE_ADDRESS_MIRROR;
    imtbl_samplers[2].Desc.AddressW        = TEXTURE_ADDRESS_MIRROR;

    imtbl_samplers[3].ShaderStages         = SHADER_TYPE_PIXEL;
    imtbl_samplers[3].SamplerOrTextureName = "g_tex2DMtrlMap";
    imtbl_samplers[3].Desc                 = imtbl_samplers[2].Desc;

    imtbl_samplers[4].ShaderStages         = SHADER_TYPE_PIXEL;
    imtbl_samplers[4].SamplerOrTextureName = "g_tex2DShadowMap";
    imtbl_samplers[4].Desc.MinFilter       = FILTER_TYPE_COMPARISON_LINEAR;
    imtbl_samplers[4].Desc.MagFilter       = FILTER_TYPE_COMPARISON_LINEAR;
    imtbl_samplers[4].Desc.MipFilter       = FILTER_TYPE_COMPARISON_LINEAR;
    imtbl_samplers[4].Desc.ComparisonFunc  = COMPARISON_FUNC_LESS;

    ShaderMacroHelper macros;
    macros.AddShaderMacro("TEXTURING_MODE", params_.texturing_mode);
    macros.AddShaderMacro("NUM_TILE_TEXTURES", NUM_TILE_TEXTURES);
    macros.AddShaderMacro("NUM_SHADOW_CASCADES", params_.num_shadow_cascades);
    macros.AddShaderMacro("BEST_CASCADE_SEARCH", params_.best_cascade_search ? true : false);
    macros.AddShaderMacro("SHADOW_FILTER_SIZE", params_.fixed_shadow_filter_size);
    macros.AddShaderMacro("FILTER_ACROSS_CASCADES", params_.filter_across_shadow_cascades);
    macros.Finalize();
    shader_create_info.Macros = macros;

    RefCntAutoPtr<IShader> hemisphere_ps;
    device_->CreateShader(shader_create_info, &hemisphere_ps);

    LayoutElement layout_inputs[] = {{0, 0, 3, VT_FLOAT32}, {1, 0, 2, VT_FLOAT32}}; //layout attribs

    GraphicsPipelineStateCreateInfo pso_create_info;
    PipelineStateDesc&              pso_desc = pso_create_info.PSODesc;

    pso_desc.Name = "RenderHemisphere";

    ShaderResourceVariableDesc shader_res_vars[] = {
      {SHADER_TYPE_VERTEX, "cbCameraAttribs", SHADER_RESOURCE_VARIABLE_TYPE_MUTABLE},
      {SHADER_TYPE_VERTEX, "cbLightAttribs", SHADER_RESOURCE_VARIABLE_TYPE_MUTABLE},
      {SHADER_TYPE_VERTEX, "cbTerrainAttribs", SHADER_RESOURCE_VARIABLE_TYPE_MUTABLE},
      {SHADER_TYPE_VERTEX, "cbParticipatingMediaScatteringParams", SHADER_RESOURCE_VARIABLE_TYPE_STATIC},
      {SHADER_TYPE_VERTEX, "g_tex2DOccludedNetDensityToAtmTop", SHADER_RESOURCE_VARIABLE_TYPE_DYNAMIC},
      {SHADER_TYPE_VERTEX, "g_tex2DAmbientSkylight", SHADER_RESOURCE_VARIABLE_TYPE_DYNAMIC},
      {SHADER_TYPE_PIXEL, "g_tex2DShadowMap", SHADER_RESOURCE_VARIABLE_TYPE_DYNAMIC},
    };
    pso_desc.ResourceLayout.Variables            = shader_res_vars;
    pso_desc.ResourceLayout.NumVariables         = _countof(shader_res_vars);
    pso_desc.ResourceLayout.ImmutableSamplers    = imtbl_samplers.data();
    pso_desc.ResourceLayout.NumImmutableSamplers = static_cast<Uint32>(imtbl_samplers.size());

    auto& graphics_pipeline = pso_create_info.GraphicsPipeline;

    graphics_pipeline.RasterizerDesc.FillMode              = FILL_MODE_SOLID;
    graphics_pipeline.RasterizerDesc.CullMode              = CULL_MODE_BACK;
    graphics_pipeline.RasterizerDesc.FrontCounterClockwise = True;
    graphics_pipeline.InputLayout.LayoutElements           = layout_inputs;
    graphics_pipeline.InputLayout.NumElements              = _countof(layout_inputs);
    pso_create_info.pVS                                     = hemisphere_vs_;
    pso_create_info.pPS                                     = hemisphere_ps;
    graphics_pipeline.RTVFormats[0]                        = params_.DstRTVFormat;
    graphics_pipeline.NumRenderTargets                     = 1;
    graphics_pipeline.DSVFormat                            = TEX_FORMAT_D32_FLOAT;
    graphics_pipeline.PrimitiveTopology                    = PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP;
    device_->CreateGraphicsPipelineState(pso_create_info, &hemisphere_pso_);
    hemisphere_pso_->BindStaticResources(SHADER_TYPE_VERTEX | SHADER_TYPE_PIXEL, resource_mapping_, BIND_SHADER_RESOURCES_VERIFY_ALL_RESOLVED);
    hemisphere_pso_->CreateShaderResourceBinding(&hemisphere_srb_, true);
    hemisphere_srb_->BindResources(SHADER_TYPE_VERTEX, resource_mapping_, BIND_SHADER_RESOURCES_KEEP_EXISTING);
  }

  ViewFrustumExt view_frustum;
  auto           dev_type = device_->GetDeviceInfo().Type;
  ExtractViewFrustumPlanesFromMatrix(camera_view_proj_matrix, view_frustum, dev_type == RENDER_DEVICE_TYPE_D3D11 || dev_type == RENDER_DEVICE_TYPE_D3D12);

  {
    MapHelper<TerrainAttribs> terrain_attribs(context, terrain_attribs_buffer_, MAP_WRITE, MAP_FLAG_DISCARD);
    *terrain_attribs = params_.terrain_attribs;
  }

  IBuffer* ppBuffers[1] = {vertex_buffer_};
  context->SetVertexBuffers(0, 1, ppBuffers, nullptr, RESOURCE_STATE_TRANSITION_MODE_TRANSITION, SET_VERTEX_BUFFERS_FLAG_RESET);

  if (z_only_pass)
  {
    context->SetPipelineState(hemisphere_z_only_pso_);
    context->CommitShaderResources(hemisphere_z_only_srb_, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);
  }
  else
  {
    shadow_map_srv->SetSampler(comparison_sampler_);
    context->SetPipelineState(hemisphere_pso_);

    hemisphere_srb_->GetVariableByName(SHADER_TYPE_PIXEL, "g_tex2DShadowMap")->Set(shadow_map_srv);
    hemisphere_srb_->GetVariableByName(SHADER_TYPE_VERTEX, "g_tex2DOccludedNetDensityToAtmTop")->Set(pre_computeD_net_density_srv);
    hemisphere_srb_->GetVariableByName(SHADER_TYPE_VERTEX, "g_tex2DAmbientSkylight")->Set(ambient_skylight_srv);

    context->CommitShaderResources(hemisphere_srb_, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);
  }

  for (auto mesh_iter = sphere_meshes_.begin(); mesh_iter != sphere_meshes_.end(); ++mesh_iter)
  {
    if (GetBoxVisibility(view_frustum, mesh_iter->bound_box, z_only_pass ? FRUSTUM_PLANE_FLAG_OPEN_NEAR : FRUSTUM_PLANE_FLAG_FULL_FRUSTUM) !=
        BoxVisibility::Invisible)
    {
      context->SetIndexBuffer(mesh_iter->index_buffer, 0, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);
      DrawIndexedAttribs draw_attrs(mesh_iter->num_indices, VT_UINT32, DRAW_FLAG_VERIFY_ALL);
      context->DrawIndexed(draw_attrs);
    }
  }
}

} // namespace Diligent
