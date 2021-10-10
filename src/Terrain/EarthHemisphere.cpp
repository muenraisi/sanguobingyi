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
      // To move from one strip to another, we have to generate two degenerate triangles by duplicating the last vertex in previous strip
      // and the first vertex in new strip
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
        indices_.push_back(base_index + index_generator_(init_col, init_row + row + static_cast<IndexType>(quad_tri_type)) + 1);
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
  StdIndexGenerator(IndexType interval) : interaval_(interval) {}

  IndexType operator()(IndexType col, IndexType row, IndexType init = 0) { return col + row * interaval_ + init; }

private:
  IndexType interaval_;
};

typedef TriStrip<Uint32, StdIndexGenerator> StdTriStrip32;

void ComputeVertexHeight(HemisphereVertex& vertex, class ElevationDataSource* elev_data_src, float sampling_step, float height_scale)
{
  float3& world_pos = vertex.world_pos;

  float hm_col     = world_pos.x / sampling_step; // sampling_step is the size by which the height map repeat
  float hm_row     = world_pos.z / sampling_step;
  float altitude = elev_data_src->GetInterpolateHeight(hm_col, hm_row);
  int   col_offset, row_offset;
  elev_data_src->GetOffsets(col_offset, row_offset);
  // +0.5f due to translate to orign£¿
  vertex.mask_uv.x = (hm_col + static_cast<float>(col_offset) + 0.5f) / static_cast<float>(elev_data_src->GetDimCol());
  vertex.mask_uv.y = (hm_row + static_cast<float>(row_offset) + 0.5f) / static_cast<float>(elev_data_src->GetDimRow());

  float3 sphere_normal = normalize(world_pos);
  world_pos += sphere_normal * altitude * height_scale;
}

class RingMeshBuilder
{
public:
  RingMeshBuilder(IRenderDevice* pDevice, const std::vector<HemisphereVertex>& VB, int iGridDimenion, std::vector<RingSectorMesh>& RingMeshes) :
    m_pDevice(pDevice), m_RingMeshes(RingMeshes), m_VB(VB), m_iGridDimenion(iGridDimenion)
  {}

  void CreateMesh(int iBaseIndex, int iStartCol, int iStartRow, int iNumCols, int iNumRows, QUAD_TRI_TYPE QuadTriangType)
  {
    m_RingMeshes.push_back(RingSectorMesh());
    auto& CurrMesh = m_RingMeshes.back();

    std::vector<Uint32> IB;
    StdTriStrip32       TriStrip(IB, StdIndexGenerator<Uint32>(m_iGridDimenion));
    TriStrip.AddStrip(iBaseIndex, iStartCol, iStartRow, iNumCols, iNumRows, QuadTriangType);

    CurrMesh.uiNumIndices = (Uint32)IB.size();

    // Prepare buffer description
    BufferDesc IndexBufferDesc;
    IndexBufferDesc.Name      = "Ring mesh index buffer";
    IndexBufferDesc.Size      = (Uint32)(IB.size() * sizeof(IB[0]));
    IndexBufferDesc.BindFlags = BIND_INDEX_BUFFER;
    IndexBufferDesc.Usage     = USAGE_IMMUTABLE;
    BufferData IBInitData;
    IBInitData.pData    = IB.data();
    IBInitData.DataSize = IndexBufferDesc.Size;
    // Create the buffer
    m_pDevice->CreateBuffer(IndexBufferDesc, &IBInitData, &CurrMesh.pIndBuff);
    VERIFY(CurrMesh.pIndBuff, "Failed to create index buffer");

    // Compute bounding box
    auto& BB = CurrMesh.BndBox;
    BB.Max   = float3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    BB.Min   = float3(+FLT_MAX, +FLT_MAX, +FLT_MAX);
    for (auto Ind = IB.begin(); Ind != IB.end(); ++Ind)
    {
      const auto& CurrVert = m_VB[*Ind].world_pos;

      BB.Min = std::min(BB.Min, CurrVert);
      BB.Max = std::max(BB.Max, CurrVert);
    }
  }

private:
  RefCntAutoPtr<IRenderDevice>         m_pDevice;
  std::vector<RingSectorMesh>&         m_RingMeshes;
  const std::vector<HemisphereVertex>& m_VB; // vertex buffer
  const size_t                         m_iGridDimenion;
};

void GenerateSphereGeometry(IRenderDevice*                 pDevice,
                            const float                    earth_radius,
                            size_t                         grid_dim, // grid dimension
                            const size_t                   num_rings,
                            class ElevationDataSource*     elevation_data_source,
                            float                          sampling_step,
                            float                          height_scale,
                            std::vector<HemisphereVertex>& vert_buff,
                            std::vector<RingSectorMesh>&   sphere_meshes)
{
  if ((grid_dim - 1) % 4 != 0)
  {
    grid_dim = RenderingParams().m_iRingDimension;
    UNEXPECTED("Grid dimension must be 4k+1. Defaulting to ",
               grid_dim); // Grid dimension of projection sphere, due to
                          // symmetry, it should be 2K+1, but why 4K+1
  }
  const size_t              grid_midst = (grid_dim - 1) / 2;
  const size_t              grid_quart = (grid_dim - 1) / 4;
  StdIndexGenerator<Uint32> std_index_generator(grid_dim);

  // const int iLargestGridScale = iGridDimension << (iNumRings-1);

  RingMeshBuilder ring_mesh_builder(pDevice, vert_buff, grid_dim, sphere_meshes);
  size_t          start_ring = 0;
  vert_buff.reserve((num_rings - start_ring) * grid_dim * grid_dim);
  for (size_t ring = start_ring; ring < num_rings; ++ring)
  {
    size_t now_grid_init = vert_buff.size();
    // Why not resize once directly?
    vert_buff.resize(vert_buff.size() + grid_dim * grid_dim);
    // What is grid scale? grid_scale = 2^-(num_rings-1-i_ring)
    float grid_scale = 1.f / (float)(1 << (num_rings - 1 - ring));
    // Fill vertex buffer
    for (size_t row = 0; row < grid_dim; ++row)
      for (size_t col = 0; col < grid_dim; ++col)
      {
        auto& now_vert = vert_buff[std_index_generator(col, row, now_grid_init)];
        auto& pos      = now_vert.world_pos;
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

        ComputeVertexHeight(now_vert, elevation_data_source, sampling_step, height_scale);
        pos.y -= earth_radius; // Translate the top to the origin
      }

    // Align vertices on the outer boundary
    if (ring < num_rings - 1)
    {
      for (size_t i = 1; i < grid_dim - 1; i += 2)
      {
        // Top & bottom boundaries
        for (size_t iRow = 0; iRow < grid_dim; iRow += grid_dim - 1)
        {
          const auto& V0 = vert_buff[now_grid_init + i - 1 + iRow * grid_dim].world_pos;
          auto&       V1 = vert_buff[now_grid_init + i + 0 + iRow * grid_dim].world_pos;
          const auto& V2 = vert_buff[now_grid_init + i + 1 + iRow * grid_dim].world_pos;
          V1             = (V0 + V2) / 2.f;
        }

        // Left & right boundaries
        for (size_t iCol = 0; iCol < grid_dim; iCol += grid_dim - 1)
        {
          const auto& V0 = vert_buff[now_grid_init + iCol + (i - 1) * grid_dim].world_pos;
          auto&       V1 = vert_buff[now_grid_init + iCol + (i + 0) * grid_dim].world_pos;
          const auto& V2 = vert_buff[now_grid_init + iCol + (i + 1) * grid_dim].world_pos;
          V1             = (V0 + V2) / 2.f;
        }
      }
    }

    // Generate indices for the current ring
    if (ring == 0)
    {
      // clang-format off
			ring_mesh_builder.CreateMesh(now_grid_init, 0, 0, grid_midst + 1, grid_midst + 1, QUAD_TRI_TYPE::DIAG_00_11);
			ring_mesh_builder.CreateMesh(now_grid_init, grid_midst, 0, grid_midst + 1, grid_midst + 1, QUAD_TRI_TYPE::DIAG_01_10);
			ring_mesh_builder.CreateMesh(now_grid_init, 0, grid_midst, grid_midst + 1, grid_midst + 1, QUAD_TRI_TYPE::DIAG_01_10);
			ring_mesh_builder.CreateMesh(now_grid_init, grid_midst, grid_midst, grid_midst + 1, grid_midst + 1, QUAD_TRI_TYPE::DIAG_00_11);
      // clang-format on
    }
    else
    {
      // clang-format off
			ring_mesh_builder.CreateMesh(now_grid_init, 0, 0, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);
			ring_mesh_builder.CreateMesh(now_grid_init, grid_quart, 0, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);

			ring_mesh_builder.CreateMesh(now_grid_init, grid_midst, 0, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);
			ring_mesh_builder.CreateMesh(now_grid_init, grid_quart * 3, 0, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);

			ring_mesh_builder.CreateMesh(now_grid_init, 0, grid_quart, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);
			ring_mesh_builder.CreateMesh(now_grid_init, 0, grid_midst, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);

			ring_mesh_builder.CreateMesh(now_grid_init, grid_quart * 3, grid_quart, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);
			ring_mesh_builder.CreateMesh(now_grid_init, grid_quart * 3, grid_midst, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);

			ring_mesh_builder.CreateMesh(now_grid_init, 0, grid_quart * 3, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);
			ring_mesh_builder.CreateMesh(now_grid_init, grid_quart, grid_quart * 3, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_01_10);

			ring_mesh_builder.CreateMesh(now_grid_init, grid_midst, grid_quart * 3, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);
			ring_mesh_builder.CreateMesh(now_grid_init, grid_quart * 3, grid_quart * 3, grid_quart + 1, grid_quart + 1, QUAD_TRI_TYPE::DIAG_00_11);
      // clang-format on
    }
  }

  // We do not need per-vertex normals as we use normal map to shade terrain
  // Sphere tangent vertex are computed in the shader
#if 0
		// Compute normals
	const float3* pV0 = nullptr;
	const float3* pV1 = &VB[IB[0]].f3WorldPos;
	const float3* pV2 = &VB[IB[1]].f3WorldPos;
	float fSign = +1;
	for (Uint32 Ind = 2; Ind < m_uiIndicesInIndBuff; ++Ind)
	{
		fSign = -fSign;
		pV0 = pV1;
		pV1 = pV2;
		pV2 = &VB[IB[Ind]].f3WorldPos;
		float3 Rib0 = *pV0 - *pV1;
		float3 Rib1 = *pV1 - *pV2;
		float3 TriN;
		D3DXVec3Cross(&TriN, &Rib0, &Rib1);
		float fLength = D3DXVec3Length(&TriN);
		if (fLength > 0.1)
		{
			TriN /= fLength * fSign;
			for (int i = -2; i <= 0; ++i)
				VB[IB[Ind + i]].f3Normal += TriN;
		}
	}
	for (auto VBIt = VB.begin(); VBIt != VB.end(); ++VBIt)
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
					auto& CurrGridN = VB[iCurrGridStart + CurrGridOffsets[Bnd] + i * iGridDimension].f3Normal;
					auto& NextGridN = VB[iNextGridStart + NextGridPffsets[Bnd] + (iGridQuart + i / 2) * iGridDimension].f3Normal;
					auto NewN = CurrGridN + NextGridN;
					D3DXVec3Normalize(&NewN, &NewN);
					CurrGridN = NextGridN = NewN;
					if (i > 1)
					{
						auto& PrevCurrGridN = VB[iCurrGridStart + CurrGridOffsets[Bnd] + (i - 2) * iGridDimension].f3Normal;
						auto MiddleN = PrevCurrGridN + NewN;
						D3DXVec3Normalize(&VB[iCurrGridStart + CurrGridOffsets[Bnd] + (i - 1) * iGridDimension].f3Normal, &MiddleN);
					}
				}

				// Bottom and top boundaries
				{
					auto& CurrGridN = VB[iCurrGridStart + i + CurrGridOffsets[Bnd] * iGridDimension].f3Normal;
					auto& NextGridN = VB[iNextGridStart + (iGridQuart + i / 2) + NextGridPffsets[Bnd] * iGridDimension].f3Normal;
					auto NewN = CurrGridN + NextGridN;
					D3DXVec3Normalize(&NewN, &NewN);
					CurrGridN = NextGridN = NewN;
					if (i > 1)
					{
						auto& PrevCurrGridN = VB[iCurrGridStart + (i - 2) + CurrGridOffsets[Bnd] * iGridDimension].f3Normal;
						auto MiddleN = PrevCurrGridN + NewN;
						D3DXVec3Normalize(&VB[iCurrGridStart + (i - 1) + CurrGridOffsets[Bnd] * iGridDimension].f3Normal, &MiddleN);
					}
				}
			}
		}
	}
#endif
}

void EarthHemsiphere::RenderNormalMap(IRenderDevice*  pDevice,
                                      IDeviceContext* pContext,
                                      const Uint16*   pHeightMap,
                                      size_t          HeightMapPitch,
                                      size_t          iHeightMapDim,
                                      ITexture*       ptex2DNormalMap)
{
  TextureDesc HeightMapDesc;
  HeightMapDesc.Name      = "Height map texture";
  HeightMapDesc.Type      = RESOURCE_DIM_TEX_2D;
  HeightMapDesc.Width     = iHeightMapDim;
  HeightMapDesc.Height    = iHeightMapDim;
  HeightMapDesc.Format    = TEX_FORMAT_R16_UINT;
  HeightMapDesc.Usage     = USAGE_IMMUTABLE;
  HeightMapDesc.BindFlags = BIND_SHADER_RESOURCE;
  HeightMapDesc.MipLevels = ComputeMipLevelsCount(HeightMapDesc.Width, HeightMapDesc.Height);

  std::vector<Uint16> CoarseMipLevels;
  CoarseMipLevels.resize(iHeightMapDim / 2 * iHeightMapDim);

  std::vector<TextureSubResData> InitData(HeightMapDesc.MipLevels);
  InitData[0].pData            = pHeightMap;
  InitData[0].Stride           = (Uint32)HeightMapPitch * sizeof(pHeightMap[0]);
  const Uint16* pFinerMipLevel = pHeightMap;
  Uint16*       pCurrMipLevel  = &CoarseMipLevels[0];
  size_t        FinerMipPitch  = HeightMapPitch;
  size_t        CurrMipPitch   = iHeightMapDim / 2;
  for (size_t uiMipLevel = 1; uiMipLevel < HeightMapDesc.MipLevels; ++uiMipLevel)
  {
    auto MipWidth  = HeightMapDesc.Width >> uiMipLevel;
    auto MipHeight = HeightMapDesc.Height >> uiMipLevel;
    for (size_t uiRow = 0; uiRow < MipHeight; ++uiRow)
    {
      for (size_t uiCol = 0; uiCol < MipWidth; ++uiCol)
      {
        int iAverageHeight = 0;
        for (size_t i = 0; i < 2; ++i)
          for (size_t j = 0; j < 2; ++j) iAverageHeight += pFinerMipLevel[(uiCol * 2 + i) + (uiRow * 2 + j) * size_t{FinerMipPitch}];
        pCurrMipLevel[uiCol + uiRow * CurrMipPitch] = (Uint16)(iAverageHeight >> 2);
      }
    }

    InitData[uiMipLevel].pData  = pCurrMipLevel;
    InitData[uiMipLevel].Stride = (Uint32)CurrMipPitch * sizeof(*pCurrMipLevel);
    pFinerMipLevel              = pCurrMipLevel;
    FinerMipPitch               = CurrMipPitch;
    pCurrMipLevel += MipHeight * CurrMipPitch;
    CurrMipPitch = iHeightMapDim / 2;
  }

  RefCntAutoPtr<ITexture> ptex2DHeightMap;
  TextureData             HeigtMapInitData;
  HeigtMapInitData.pSubResources   = InitData.data();
  HeigtMapInitData.NumSubresources = (Uint32)InitData.size();
  pDevice->CreateTexture(HeightMapDesc, &HeigtMapInitData, &ptex2DHeightMap);
  VERIFY(ptex2DHeightMap, "Failed to create height map texture");

  m_pResMapping->AddResource("g_tex2DElevationMap", ptex2DHeightMap->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE), true);

  RefCntAutoPtr<IBuffer> pcbNMGenerationAttribs;
  CreateUniformBuffer(pDevice, sizeof(NMGenerationAttribs), "NM Generation Attribs CB", &pcbNMGenerationAttribs);

  m_pResMapping->AddResource("cbNMGenerationAttribs", pcbNMGenerationAttribs, true);

  RefCntAutoPtr<IShaderSourceInputStreamFactory> pShaderSourceFactory;
  m_pDevice->GetEngineFactory()->CreateDefaultShaderSourceStreamFactory("shaders\\;shaders\\terrain", &pShaderSourceFactory);

  ShaderCreateInfo ShaderCI;
  ShaderCI.pShaderSourceStreamFactory = pShaderSourceFactory;
  ShaderCI.FilePath                   = "ScreenSizeQuadVS.fx";
  ShaderCI.EntryPoint                 = "GenerateScreenSizeQuadVS";
  ShaderCI.SourceLanguage             = SHADER_SOURCE_LANGUAGE_HLSL;
  ShaderCI.UseCombinedTextureSamplers = true;
  ShaderCI.Desc.ShaderType            = SHADER_TYPE_VERTEX;
  ShaderCI.Desc.Name                  = "GenerateScreenSizeQuadVS";
  RefCntAutoPtr<IShader> pScreenSizeQuadVS;
  pDevice->CreateShader(ShaderCI, &pScreenSizeQuadVS);

  ShaderCI.FilePath        = "GenerateNormalMapPS.fx";
  ShaderCI.EntryPoint      = "GenerateNormalMapPS";
  ShaderCI.Desc.ShaderType = SHADER_TYPE_PIXEL;
  ShaderCI.Desc.Name       = "GenerateNormalMapPS";

  RefCntAutoPtr<IShader> pGenerateNormalMapPS;
  pDevice->CreateShader(ShaderCI, &pGenerateNormalMapPS);

  GraphicsPipelineStateCreateInfo PSOCreateInfo;
  PipelineStateDesc&              PSODesc = PSOCreateInfo.PSODesc;

  PSODesc.ResourceLayout.DefaultVariableType = SHADER_RESOURCE_VARIABLE_TYPE_STATIC;
  // clang-format off
    ShaderResourceVariableDesc ShaderVars[] = 
    {
        {SHADER_TYPE_PIXEL, "g_tex2DElevationMap",   SHADER_RESOURCE_VARIABLE_TYPE_STATIC},
        {SHADER_TYPE_PIXEL, "cbNMGenerationAttribs", SHADER_RESOURCE_VARIABLE_TYPE_STATIC}
    };
  // clang-format on

  PSODesc.ResourceLayout.NumVariables = _countof(ShaderVars);
  PSODesc.ResourceLayout.Variables    = ShaderVars;

  PSODesc.Name           = "Render Normal Map";
  auto& GraphicsPipeline = PSOCreateInfo.GraphicsPipeline;

  GraphicsPipeline.DepthStencilDesc.DepthEnable         = false;
  GraphicsPipeline.DepthStencilDesc.DepthWriteEnable    = false;
  GraphicsPipeline.RasterizerDesc.FillMode              = FILL_MODE_SOLID;
  GraphicsPipeline.RasterizerDesc.CullMode              = CULL_MODE_NONE;
  GraphicsPipeline.RasterizerDesc.FrontCounterClockwise = true;
  PSOCreateInfo.pVS                                     = pScreenSizeQuadVS;
  PSOCreateInfo.pPS                                     = pGenerateNormalMapPS;
  GraphicsPipeline.NumRenderTargets                     = 1;
  GraphicsPipeline.RTVFormats[0]                        = TEX_FORMAT_RG8_UNORM;
  GraphicsPipeline.PrimitiveTopology                    = PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP;
  RefCntAutoPtr<IPipelineState> pRenderNormalMapPSO;
  pDevice->CreateGraphicsPipelineState(PSOCreateInfo, &pRenderNormalMapPSO);
  pRenderNormalMapPSO->BindStaticResources(SHADER_TYPE_VERTEX | SHADER_TYPE_PIXEL, m_pResMapping, BIND_SHADER_RESOURCES_VERIFY_ALL_RESOLVED);

  RefCntAutoPtr<IShaderResourceBinding> pRenderNormalMapSRB;
  pRenderNormalMapPSO->CreateShaderResourceBinding(&pRenderNormalMapSRB, true);

  pContext->SetPipelineState(pRenderNormalMapPSO);
  pContext->CommitShaderResources(pRenderNormalMapSRB, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);

  const auto& NormalMapDesc = ptex2DNormalMap->GetDesc();
  for (Uint32 uiMipLevel = 0; uiMipLevel < NormalMapDesc.MipLevels; ++uiMipLevel)
  {
    TextureViewDesc TexViewDesc;
    TexViewDesc.ViewType        = TEXTURE_VIEW_RENDER_TARGET;
    TexViewDesc.MostDetailedMip = uiMipLevel;
    RefCntAutoPtr<ITextureView> ptex2DNormalMapRTV;
    ptex2DNormalMap->CreateView(TexViewDesc, &ptex2DNormalMapRTV);

    ITextureView* pRTVs[] = {ptex2DNormalMapRTV};
    pContext->SetRenderTargets(_countof(pRTVs), pRTVs, nullptr, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);

    {
      MapHelper<NMGenerationAttribs> NMGenerationAttribs(pContext, pcbNMGenerationAttribs, MAP_WRITE, MAP_FLAG_DISCARD);
      NMGenerationAttribs->m_fElevationScale        = m_Params.m_TerrainAttribs.m_fElevationScale;
      NMGenerationAttribs->m_fSampleSpacingInterval = m_Params.m_TerrainAttribs.m_fElevationSamplingInterval;
      NMGenerationAttribs->m_iMIPLevel              = static_cast<int>(uiMipLevel);
    }

    DrawAttribs DrawAttrs(4, DRAW_FLAG_VERIFY_ALL);
    pContext->Draw(DrawAttrs);
  }

  // Remove elevation map from resource mapping to release the resource
  m_pResMapping->RemoveResourceByName("g_tex2DElevationMap");
}

void EarthHemsiphere::Create(class ElevationDataSource* pDataSource,
                             const RenderingParams&     Params,
                             IRenderDevice*             pDevice,
                             IDeviceContext*            pContext,
                             const Char*                MaterialMaskPath,
                             const Char*                TileTexturePath[],
                             const Char*                TileNormalMapPath[],
                             IBuffer*                   pcbCameraAttribs,
                             IBuffer*                   pcbLightAttribs,
                             IBuffer*                   pcMediaScatteringParams)
{
  m_Params  = Params;
  m_pDevice = pDevice;

  const Uint16* pHeightMap;
  size_t        HeightMapPitch;
  pDataSource->GetDataPtr(pHeightMap, HeightMapPitch);
  Uint32 iHeightMapDim = pDataSource->GetDimCol();
  VERIFY_EXPR(iHeightMapDim == pDataSource->GetDimRow());

  TextureDesc NormalMapDesc;
  NormalMapDesc.Name      = "Normal map texture";
  NormalMapDesc.Type      = RESOURCE_DIM_TEX_2D;
  NormalMapDesc.Width     = iHeightMapDim;
  NormalMapDesc.Height    = iHeightMapDim;
  NormalMapDesc.Format    = TEX_FORMAT_RG8_UNORM;
  NormalMapDesc.Usage     = USAGE_DEFAULT;
  NormalMapDesc.BindFlags = BIND_SHADER_RESOURCE | BIND_RENDER_TARGET;
  NormalMapDesc.MipLevels = 0;

  RefCntAutoPtr<ITexture> ptex2DNormalMap;
  pDevice->CreateTexture(NormalMapDesc, nullptr, &ptex2DNormalMap);
  m_ptex2DNormalMapSRV = ptex2DNormalMap->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE);

  CreateUniformBuffer(pDevice, sizeof(TerrainAttribs), "Terrain Attribs CB", &m_pcbTerrainAttribs);

  ResourceMappingDesc ResMappingDesc;
  // clang-format off
    ResourceMappingEntry pEntries[] = 
    { 
        { "cbCameraAttribs", pcbCameraAttribs }, 
        { "cbTerrainAttribs", m_pcbTerrainAttribs}, 
        { "cbLightAttribs", pcbLightAttribs}, 
        { "g_tex2DNormalMap", ptex2DNormalMap->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE) }, 
        { "cbParticipatingMediaScatteringParams", pcMediaScatteringParams },
        {} 
    };
  // clang-format on
  ResMappingDesc.pEntries = pEntries;
  pDevice->CreateResourceMapping(ResMappingDesc, &m_pResMapping);

  RefCntAutoPtr<ITexture> ptex2DMtrlMask;
  CreateTextureFromFile(MaterialMaskPath, TextureLoadInfo(), pDevice, &ptex2DMtrlMask);
  auto ptex2DMtrlMaskSRV = ptex2DMtrlMask->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE);
  m_pResMapping->AddResource("g_tex2DMtrlMap", ptex2DMtrlMaskSRV, true);

  // Load tiles
  IDeviceObject*          ptex2DTileDiffuseSRV[NUM_TILE_TEXTURES] = {};
  RefCntAutoPtr<ITexture> ptex2DTileDiffuse[NUM_TILE_TEXTURES];
  IDeviceObject*          ptex2DTileNMSRV[NUM_TILE_TEXTURES] = {};
  RefCntAutoPtr<ITexture> ptex2DTileNM[NUM_TILE_TEXTURES];
  for (size_t iTileTex = 0; iTileTex < static_cast<size_t>(NUM_TILE_TEXTURES); iTileTex++)
  {
    {
      TextureLoadInfo DiffMapLoadInfo;
      DiffMapLoadInfo.IsSRGB = false;
      CreateTextureFromFile(TileTexturePath[iTileTex], DiffMapLoadInfo, pDevice, &ptex2DTileDiffuse[iTileTex]);
      ptex2DTileDiffuseSRV[iTileTex] = ptex2DTileDiffuse[iTileTex]->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE);
    }

    {
      CreateTextureFromFile(TileNormalMapPath[iTileTex], TextureLoadInfo(), pDevice, &ptex2DTileNM[iTileTex]);
      ptex2DTileNMSRV[iTileTex] = ptex2DTileNM[iTileTex]->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE);
    }
  }
  m_pResMapping->AddResourceArray("g_tex2DTileDiffuse", 0, ptex2DTileDiffuseSRV, NUM_TILE_TEXTURES, true);
  m_pResMapping->AddResourceArray("g_tex2DTileNM", 0, ptex2DTileNMSRV, NUM_TILE_TEXTURES, true);

  m_pDevice->CreateSampler(Sam_ComparsionLinearClamp, &m_pComparisonSampler);

  RenderNormalMap(pDevice, pContext, pHeightMap, HeightMapPitch, iHeightMapDim, ptex2DNormalMap);

  RefCntAutoPtr<IShaderSourceInputStreamFactory> pShaderSourceFactory;
  m_pDevice->GetEngineFactory()->CreateDefaultShaderSourceStreamFactory("shaders;shaders\\terrain;", &pShaderSourceFactory);

  {
    ShaderCreateInfo ShaderCI;
    ShaderCI.pShaderSourceStreamFactory = pShaderSourceFactory;
    ShaderCI.FilePath                   = "HemisphereVS.fx";
    ShaderCI.EntryPoint                 = "HemisphereVS";
    ShaderCI.SourceLanguage             = SHADER_SOURCE_LANGUAGE_HLSL;
    ShaderCI.UseCombinedTextureSamplers = true;
    ShaderCI.Desc.ShaderType            = SHADER_TYPE_VERTEX;
    ShaderCI.Desc.Name                  = "HemisphereVS";

    pDevice->CreateShader(ShaderCI, &m_pHemisphereVS);
  }

  {
    ShaderCreateInfo ShaderCI;
    ShaderCI.pShaderSourceStreamFactory = pShaderSourceFactory;
    ShaderCI.FilePath                   = "HemisphereZOnlyVS.fx";
    ShaderCI.EntryPoint                 = "HemisphereZOnlyVS";
    ShaderCI.SourceLanguage             = SHADER_SOURCE_LANGUAGE_HLSL;
    ShaderCI.UseCombinedTextureSamplers = true;
    ShaderCI.Desc.ShaderType            = SHADER_TYPE_VERTEX;
    ShaderCI.Desc.Name                  = "HemisphereZOnlyVS";
    RefCntAutoPtr<IShader> pHemisphereZOnlyVS;
    pDevice->CreateShader(ShaderCI, &pHemisphereZOnlyVS);

    GraphicsPipelineStateCreateInfo PSOCreateInfo;
    PipelineStateDesc&              PSODesc = PSOCreateInfo.PSODesc;

    PSODesc.Name                                          = "Render Hemisphere Z Only";
    PSODesc.ResourceLayout.DefaultVariableType            = SHADER_RESOURCE_VARIABLE_TYPE_STATIC;
    auto& GraphicsPipeline                                = PSOCreateInfo.GraphicsPipeline;
    GraphicsPipeline.DepthStencilDesc                     = DSS_Default;
    GraphicsPipeline.RasterizerDesc.FillMode              = FILL_MODE_SOLID;
    GraphicsPipeline.RasterizerDesc.CullMode              = CULL_MODE_BACK;
    GraphicsPipeline.RasterizerDesc.DepthClipEnable       = False;
    GraphicsPipeline.RasterizerDesc.FrontCounterClockwise = True;
    // clang-format off
        LayoutElement Inputs[] =
        {
            {0, 0, 3, VT_FLOAT32, False, 0, (3+2)*4}
        };
    // clang-format on
    GraphicsPipeline.InputLayout.LayoutElements = Inputs;
    GraphicsPipeline.InputLayout.NumElements    = _countof(Inputs);
    GraphicsPipeline.DSVFormat                  = m_Params.ShadowMapFormat;
    GraphicsPipeline.PrimitiveTopology          = PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP;
    PSOCreateInfo.pVS                           = pHemisphereZOnlyVS;
    pDevice->CreateGraphicsPipelineState(PSOCreateInfo, &m_pHemisphereZOnlyPSO);
    m_pHemisphereZOnlyPSO->BindStaticResources(SHADER_TYPE_VERTEX | SHADER_TYPE_PIXEL, m_pResMapping, BIND_SHADER_RESOURCES_VERIFY_ALL_RESOLVED);
    m_pHemisphereZOnlyPSO->CreateShaderResourceBinding(&m_pHemisphereZOnlySRB, true);
  }

  std::vector<HemisphereVertex> VB;
  GenerateSphereGeometry(pDevice, Diligent::AirScatteringAttribs().fEarthRadius, m_Params.m_iRingDimension, m_Params.m_iNumRings, pDataSource,
                         m_Params.m_TerrainAttribs.m_fElevationSamplingInterval, m_Params.m_TerrainAttribs.m_fElevationScale, VB, m_SphereMeshes);

  BufferDesc VBDesc;
  VBDesc.Name      = "Hemisphere vertex buffer";
  VBDesc.Size      = static_cast<Uint64>(VB.size() * sizeof(VB[0]));
  VBDesc.Usage     = USAGE_IMMUTABLE;
  VBDesc.BindFlags = BIND_VERTEX_BUFFER;
  BufferData VBInitData;
  VBInitData.pData    = VB.data();
  VBInitData.DataSize = VBDesc.Size;
  pDevice->CreateBuffer(VBDesc, &VBInitData, &m_pVertBuff);
  VERIFY(m_pVertBuff, "Failed to create VB");
}

void EarthHemsiphere::Render(IDeviceContext*        pContext,
                             const RenderingParams& NewParams,
                             const float3&          vCameraPosition,
                             const float4x4&        CameraViewProjMatrix,
                             ITextureView*          pShadowMapSRV,
                             ITextureView*          pPrecomputedNetDensitySRV,
                             ITextureView*          pAmbientSkylightSRV,
                             bool                   bZOnlyPass)
{
  if (m_Params.m_iNumShadowCascades != NewParams.m_iNumShadowCascades || m_Params.m_bBestCascadeSearch != NewParams.m_bBestCascadeSearch ||
      m_Params.m_FilterAcrossShadowCascades != NewParams.m_FilterAcrossShadowCascades ||
      m_Params.m_FixedShadowFilterSize != NewParams.m_FixedShadowFilterSize || m_Params.DstRTVFormat != NewParams.DstRTVFormat)
  {
    m_pHemispherePSO.Release();
    m_pHemisphereSRB.Release();
  }

  m_Params = NewParams;

#if 0
    if( GetAsyncKeyState(VK_F9) )
    {
        m_RenderEarthHemisphereTech.Release();
    }
#endif

  if (!m_pHemispherePSO)
  {
    ShaderCreateInfo Attrs;
    Attrs.FilePath                   = "HemispherePS.fx";
    Attrs.EntryPoint                 = "HemispherePS";
    Attrs.Desc.ShaderType            = SHADER_TYPE_PIXEL;
    Attrs.Desc.Name                  = "HemispherePS";
    Attrs.UseCombinedTextureSamplers = true;
    Attrs.SourceLanguage             = SHADER_SOURCE_LANGUAGE_HLSL;
    RefCntAutoPtr<IShaderSourceInputStreamFactory> pShaderSourceFactory;
    m_pDevice->GetEngineFactory()->CreateDefaultShaderSourceStreamFactory("shaders;shaders\\terrain;", &pShaderSourceFactory);
    Attrs.pShaderSourceStreamFactory = pShaderSourceFactory;

    std::array<ImmutableSamplerDesc, 5> ImtblSamplers = {};

    ImtblSamplers[0].ShaderStages         = SHADER_TYPE_PIXEL;
    ImtblSamplers[0].SamplerOrTextureName = "g_tex2DTileDiffuse";
    ImtblSamplers[0].Desc.AddressU        = TEXTURE_ADDRESS_WRAP;
    ImtblSamplers[0].Desc.AddressV        = TEXTURE_ADDRESS_WRAP;
    ImtblSamplers[0].Desc.AddressW        = TEXTURE_ADDRESS_WRAP;

    ImtblSamplers[1].ShaderStages         = SHADER_TYPE_PIXEL;
    ImtblSamplers[1].SamplerOrTextureName = "g_tex2DTileNM";
    ImtblSamplers[1].Desc                 = ImtblSamplers[0].Desc;

    ImtblSamplers[2].ShaderStages         = SHADER_TYPE_PIXEL;
    ImtblSamplers[2].SamplerOrTextureName = "g_tex2DNormalMap";
    ImtblSamplers[2].Desc.AddressU        = TEXTURE_ADDRESS_MIRROR;
    ImtblSamplers[2].Desc.AddressV        = TEXTURE_ADDRESS_MIRROR;
    ImtblSamplers[2].Desc.AddressW        = TEXTURE_ADDRESS_MIRROR;

    ImtblSamplers[3].ShaderStages         = SHADER_TYPE_PIXEL;
    ImtblSamplers[3].SamplerOrTextureName = "g_tex2DMtrlMap";
    ImtblSamplers[3].Desc                 = ImtblSamplers[2].Desc;

    ImtblSamplers[4].ShaderStages         = SHADER_TYPE_PIXEL;
    ImtblSamplers[4].SamplerOrTextureName = "g_tex2DShadowMap";
    ImtblSamplers[4].Desc.MinFilter       = FILTER_TYPE_COMPARISON_LINEAR;
    ImtblSamplers[4].Desc.MagFilter       = FILTER_TYPE_COMPARISON_LINEAR;
    ImtblSamplers[4].Desc.MipFilter       = FILTER_TYPE_COMPARISON_LINEAR;
    ImtblSamplers[4].Desc.ComparisonFunc  = COMPARISON_FUNC_LESS;

    ShaderMacroHelper Macros;
    Macros.AddShaderMacro("TEXTURING_MODE", m_Params.m_TexturingMode);
    Macros.AddShaderMacro("NUM_TILE_TEXTURES", NUM_TILE_TEXTURES);
    Macros.AddShaderMacro("NUM_SHADOW_CASCADES", m_Params.m_iNumShadowCascades);
    Macros.AddShaderMacro("BEST_CASCADE_SEARCH", m_Params.m_bBestCascadeSearch ? true : false);
    Macros.AddShaderMacro("SHADOW_FILTER_SIZE", m_Params.m_FixedShadowFilterSize);
    Macros.AddShaderMacro("FILTER_ACROSS_CASCADES", m_Params.m_FilterAcrossShadowCascades);
    Macros.Finalize();
    Attrs.Macros = Macros;

    RefCntAutoPtr<IShader> pHemispherePS;
    m_pDevice->CreateShader(Attrs, &pHemispherePS);

    LayoutElement Inputs[] = {{0, 0, 3, VT_FLOAT32}, {1, 0, 2, VT_FLOAT32}};

    GraphicsPipelineStateCreateInfo PSOCreateInfo;
    PipelineStateDesc&              PSODesc = PSOCreateInfo.PSODesc;

    PSODesc.Name = "RenderHemisphere";

    ShaderResourceVariableDesc Vars[] = {
      {SHADER_TYPE_VERTEX, "cbCameraAttribs", SHADER_RESOURCE_VARIABLE_TYPE_MUTABLE},
      {SHADER_TYPE_VERTEX, "cbLightAttribs", SHADER_RESOURCE_VARIABLE_TYPE_MUTABLE},
      {SHADER_TYPE_VERTEX, "cbTerrainAttribs", SHADER_RESOURCE_VARIABLE_TYPE_MUTABLE},
      {SHADER_TYPE_VERTEX, "cbParticipatingMediaScatteringParams", SHADER_RESOURCE_VARIABLE_TYPE_STATIC},
      {SHADER_TYPE_VERTEX, "g_tex2DOccludedNetDensityToAtmTop", SHADER_RESOURCE_VARIABLE_TYPE_DYNAMIC},
      {SHADER_TYPE_VERTEX, "g_tex2DAmbientSkylight", SHADER_RESOURCE_VARIABLE_TYPE_DYNAMIC},
      {SHADER_TYPE_PIXEL, "g_tex2DShadowMap", SHADER_RESOURCE_VARIABLE_TYPE_DYNAMIC},
    };
    PSODesc.ResourceLayout.Variables            = Vars;
    PSODesc.ResourceLayout.NumVariables         = _countof(Vars);
    PSODesc.ResourceLayout.ImmutableSamplers    = ImtblSamplers.data();
    PSODesc.ResourceLayout.NumImmutableSamplers = static_cast<Uint32>(ImtblSamplers.size());

    auto& GraphicsPipeline = PSOCreateInfo.GraphicsPipeline;

    GraphicsPipeline.RasterizerDesc.FillMode              = FILL_MODE_SOLID;
    GraphicsPipeline.RasterizerDesc.CullMode              = CULL_MODE_BACK;
    GraphicsPipeline.RasterizerDesc.FrontCounterClockwise = True;
    GraphicsPipeline.InputLayout.LayoutElements           = Inputs;
    GraphicsPipeline.InputLayout.NumElements              = _countof(Inputs);
    PSOCreateInfo.pVS                                     = m_pHemisphereVS;
    PSOCreateInfo.pPS                                     = pHemispherePS;
    GraphicsPipeline.RTVFormats[0]                        = m_Params.DstRTVFormat;
    GraphicsPipeline.NumRenderTargets                     = 1;
    GraphicsPipeline.DSVFormat                            = TEX_FORMAT_D32_FLOAT;
    GraphicsPipeline.PrimitiveTopology                    = PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP;
    m_pDevice->CreateGraphicsPipelineState(PSOCreateInfo, &m_pHemispherePSO);
    m_pHemispherePSO->BindStaticResources(SHADER_TYPE_VERTEX | SHADER_TYPE_PIXEL, m_pResMapping, BIND_SHADER_RESOURCES_VERIFY_ALL_RESOLVED);
    m_pHemispherePSO->CreateShaderResourceBinding(&m_pHemisphereSRB, true);
    m_pHemisphereSRB->BindResources(SHADER_TYPE_VERTEX, m_pResMapping, BIND_SHADER_RESOURCES_KEEP_EXISTING);
  }

  ViewFrustumExt ViewFrustum;
  auto           DevType = m_pDevice->GetDeviceInfo().Type;
  ExtractViewFrustumPlanesFromMatrix(CameraViewProjMatrix, ViewFrustum, DevType == RENDER_DEVICE_TYPE_D3D11 || DevType == RENDER_DEVICE_TYPE_D3D12);

  {
    MapHelper<TerrainAttribs> TerrainAttribs(pContext, m_pcbTerrainAttribs, MAP_WRITE, MAP_FLAG_DISCARD);
    *TerrainAttribs = m_Params.m_TerrainAttribs;
  }

#if 0
    ID3D11ShaderResourceView *pSRVs[3 + 2*NUM_TILE_TEXTURES] = 
    {
        m_ptex2DNormalMapSRV,
        m_ptex2DMtrlMaskSRV,
        pShadowMapSRV
    };
    for(int iTileTex = 0; iTileTex < NUM_TILE_TEXTURES; iTileTex++)
    {
        pSRVs[3+iTileTex] = m_ptex2DTilesSRV[iTileTex];
        pSRVs[3+NUM_TILE_TEXTURES+iTileTex] = m_ptex2DTilNormalMapsSRV[iTileTex];
    }
    pd3dImmediateContext->PSSetShaderResources(1, _countof(pSRVs), pSRVs);
    pSRVs[0] = pPrecomputedNetDensitySRV;
    pSRVs[1] = pAmbientSkylightSRV;
    pd3dImmediateContext->VSSetShaderResources(0, 2, pSRVs);

    ID3D11SamplerState *pSamplers[] = {m_psamLinearMirror, m_psamLinearWrap, m_psamComaprison, m_psamLinearClamp};
	pd3dImmediateContext->VSSetSamplers(0, _countof(pSamplers), pSamplers);
	pd3dImmediateContext->PSSetSamplers(0, _countof(pSamplers), pSamplers);
#endif

  IBuffer* ppBuffers[1] = {m_pVertBuff};
  pContext->SetVertexBuffers(0, 1, ppBuffers, nullptr, RESOURCE_STATE_TRANSITION_MODE_TRANSITION, SET_VERTEX_BUFFERS_FLAG_RESET);

  if (bZOnlyPass)
  {
    pContext->SetPipelineState(m_pHemisphereZOnlyPSO);
    pContext->CommitShaderResources(m_pHemisphereZOnlySRB, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);
  }
  else
  {
    pShadowMapSRV->SetSampler(m_pComparisonSampler);
    pContext->SetPipelineState(m_pHemispherePSO);

    m_pHemisphereSRB->GetVariableByName(SHADER_TYPE_VERTEX, "g_tex2DOccludedNetDensityToAtmTop")->Set(pPrecomputedNetDensitySRV);
    m_pHemisphereSRB->GetVariableByName(SHADER_TYPE_VERTEX, "g_tex2DAmbientSkylight")->Set(pAmbientSkylightSRV);
    m_pHemisphereSRB->GetVariableByName(SHADER_TYPE_PIXEL, "g_tex2DShadowMap")->Set(pShadowMapSRV);

    pContext->CommitShaderResources(m_pHemisphereSRB, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);
  }

  for (auto MeshIt = m_SphereMeshes.begin(); MeshIt != m_SphereMeshes.end(); ++MeshIt)
  {
    if (GetBoxVisibility(ViewFrustum, MeshIt->BndBox, bZOnlyPass ? FRUSTUM_PLANE_FLAG_OPEN_NEAR : FRUSTUM_PLANE_FLAG_FULL_FRUSTUM) !=
        BoxVisibility::Invisible)
    {
      pContext->SetIndexBuffer(MeshIt->pIndBuff, 0, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);
      DrawIndexedAttribs DrawAttrs(MeshIt->uiNumIndices, VT_UINT32, DRAW_FLAG_VERIFY_ALL);
      pContext->DrawIndexed(DrawAttrs);
    }
  }
}

} // namespace Diligent
