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

#pragma once

#include "SampleBase.hpp"
#include "BasicMath.hpp"
#include "EarthHemisphere.hpp"
#include "ElevationDataSource.hpp"
#include "EpipolarLightScattering.hpp"
#include "ShadowMapManager.hpp"

namespace Diligent
{

class AtmosphereSample final : public SampleBase
{
public:
  AtmosphereSample();
  ~AtmosphereSample();

  virtual void ModifyEngineInitInfo(const ModifyEngineInitInfoAttribs& attribs) override final;

  virtual void Initialize(const SampleInitInfo& info) override final;
  virtual void Render() override final;
  virtual void Update(double CurrTime, double ElapsedTime) override final;
  virtual void WindowResize(Uint32 Width, Uint32 Height) override final;

  virtual const Char* GetSampleName() const override final { return "±øÞÄ¡¤Èý¹ú"; }

private:
  void UpdateUI();
  void CreateShadowMap();
  void RenderShadowMap(IDeviceContext* pContext, LightAttribs& LightAttribs, const float4x4& mCameraView,
                       const float4x4& mCameraProj);

  float3 m_f3LightDir = {-0.554699242f, -0.0599640049f, -0.829887390f};

  Quaternion m_CameraRotation = {0, 0, 0, 1};
  float3     m_f3CameraPos    = {0, 8000.f, 0};
  float4x4   m_mCameraView;
  float4x4   m_mCameraProj;

  RefCntAutoPtr<IBuffer> camera_attribs_;
  RefCntAutoPtr<IBuffer> light_attribs_;

  ShadowMapManager shadow_mag_manager_;
  struct ShadowSettings
  {
    Uint32 Resolution                 = 1024;
    float  fCascadePartitioningFactor = 0.95f;
    bool   bVisualizeCascades         = false;
    int    iFixedFilterSize           = 5;
  } shadow_settings_;

  RefCntAutoPtr<ISampler> comparison_sampler_;

  RenderingParams                terrain_render_params_;
  EpipolarLightScatteringAttribs epipolar_light_scattering_attribs_;

  String height_map_path_;
  String mtrl_mask_path_;
  String tile_diffuse_tex_paths_[EarthHemsiphere::NUM_TILE_TEXTURES];
  String tile_normal_tex_paths_[EarthHemsiphere::NUM_TILE_TEXTURES];

  float min_elevation_ = 0, max_elevation_ = 0;

  std::unique_ptr<ElevationDataSource> elev_data_source_;
  EarthHemsiphere                      earth_hemisphere_;
  bool                                 is_gl_device_ = false;

  std::unique_ptr<EpipolarLightScattering> epipolar_light_scattering_;

  bool   m_bEnableLightScattering = true;
  float  m_fElapsedTime           = 0.f;
  float3 custom_rlgh_beta_, custom_mie_beta_, custom_ozone_absorption_;

  RefCntAutoPtr<ITexture> m_pOffscreenColorBuffer;
  RefCntAutoPtr<ITexture> m_pOffscreenDepthBuffer;

  float      m_fCameraYaw   = 0.23f;
  float      m_fCameraPitch = 0.18f;
  MouseState m_LastMouseState;

  bool rg16u_fmt_supported = false;
  bool rg32f_fmt_supported = false;
};

} // namespace Diligent
