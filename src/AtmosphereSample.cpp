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

#include <algorithm>
#include <array>
#include <cmath>

#include "AtmosphereSample.hpp"
#include "GraphicsUtilities.h"
#include "ImGuiUtils.hpp"
#include "MapHelper.hpp"
#include "PlatformMisc.hpp"
#include "imGuIZMO.h"
#include "imgui.h"

namespace Diligent
{

SampleBase* CreateSample() { return new AtmosphereSample(); }

AtmosphereSample::AtmosphereSample() {}

void AtmosphereSample::ModifyEngineInitInfo(const ModifyEngineInitInfoAttribs& attribs)
{
  SampleBase::ModifyEngineInitInfo(attribs);

  attribs.EngineCI.Features.ComputeShaders = DEVICE_FEATURE_STATE_ENABLED;
  attribs.EngineCI.Features.DepthClamp     = DEVICE_FEATURE_STATE_OPTIONAL;
}

void AtmosphereSample::Initialize(const SampleInitInfo& info)
{
  SampleBase::Initialize(info);

  const auto& deviceInfo  = info.pDevice->GetDeviceInfo();
  is_gl_device_           = deviceInfo.IsGLDevice();
  const auto adapter_type = info.pDevice->GetAdapterInfo().Type;
  if (adapter_type == ADAPTER_TYPE_INTEGRATED)
  {
    shadow_settings_.Resolution                                = 512;
    terrain_render_params_.filter_across_shadow_cascades       = false;
    shadow_settings_.iFixedFilterSize                          = 3;
    epipolar_light_scattering_attribs_.iFirstCascadeToRayMarch = 2;
    epipolar_light_scattering_attribs_.iSingleScatteringMode   = SINGLE_SCTR_MODE_LUT;
    terrain_render_params_.num_shadow_cascades                 = 4;
    terrain_render_params_.num_rings                           = 10;
    terrain_render_params_.texturing_mode                      = RenderingParams::TM_MATERIAL_MASK;
  }

  const auto& rg16u_attribs = m_pDevice->GetTextureFormatInfoExt(TEX_FORMAT_RG16_UNORM);
  const auto& rg32f_attribs = m_pDevice->GetTextureFormatInfoExt(TEX_FORMAT_RG32_FLOAT);
  rg16u_fmt_supported     = rg16u_attribs.Supported && (rg16u_attribs.BindFlags & BIND_RENDER_TARGET);
  rg32f_fmt_supported     = rg32f_attribs.Supported && (rg32f_attribs.BindFlags & BIND_RENDER_TARGET);
  if (!rg16u_fmt_supported && !rg32f_fmt_supported) { epipolar_light_scattering_attribs_.bUse1DMinMaxTree = FALSE; }
  else
  {
    if (rg16u_fmt_supported && !rg32f_fmt_supported) epipolar_light_scattering_attribs_.bIs32BitMinMaxMipMap = FALSE;

    if (!rg16u_fmt_supported && rg32f_fmt_supported) epipolar_light_scattering_attribs_.bIs32BitMinMaxMipMap = TRUE;
  }

  custom_rlgh_beta_        = epipolar_light_scattering_attribs_.f4CustomRlghBeta;
  custom_mie_beta_         = epipolar_light_scattering_attribs_.f4CustomMieBeta;
  custom_ozone_absorption_ = epipolar_light_scattering_attribs_.f4CustomOzoneAbsorption;

  height_map_path_           = "Terrain\\HeightMap.tif";
  mtrl_mask_path_            = "Terrain\\Mask.png";
  tile_diffuse_tex_paths_[0] = "Terrain\\Tiles\\gravel_DM.dds";
  tile_diffuse_tex_paths_[1] = "Terrain\\Tiles\\grass_DM.dds";
  tile_diffuse_tex_paths_[2] = "Terrain\\Tiles\\cliff_DM.dds";
  tile_diffuse_tex_paths_[3] = "Terrain\\Tiles\\snow_DM.dds";
  tile_diffuse_tex_paths_[4] = "Terrain\\Tiles\\grassDark_DM.dds";
  tile_normal_tex_paths_[0]  = "Terrain\\Tiles\\gravel_NM.dds";
  tile_normal_tex_paths_[1]  = "Terrain\\Tiles\\grass_NM.dds";
  tile_normal_tex_paths_[2]  = "Terrain\\Tiles\\cliff_NM.dds";
  tile_normal_tex_paths_[3]  = "Terrain\\Tiles\\Snow_NM.jpg";
  tile_normal_tex_paths_[4]  = "Terrain\\Tiles\\grass_NM.dds";

  // Create data source
  try
  {
    elev_data_source_.reset(new ElevationDataSource(height_map_path_.c_str()));
    elev_data_source_->SetOffsets(terrain_render_params_.m_iColOffset, terrain_render_params_.m_iRowOffset);
    min_elevation_ = elev_data_source_->GetGlobalMinElevation() * terrain_render_params_.terrain_attribs.height_scale;
    max_elevation_ = elev_data_source_->GetGlobalMaxElevation() * terrain_render_params_.terrain_attribs.height_scale;
  }
  catch (const std::exception&)
  {
    LOG_ERROR("Failed to create elevation data source");
    return;
  }

  const Char *tile_diffuse_tex_paths[EarthHemsiphere::NUM_TILE_TEXTURES], *tile_normal_tex_paths[EarthHemsiphere::NUM_TILE_TEXTURES];
  for (size_t iTile = 0; iTile < _countof(tile_diffuse_tex_paths); ++iTile)
  {
    tile_diffuse_tex_paths[iTile]   = tile_diffuse_tex_paths_[iTile].c_str();
    tile_normal_tex_paths[iTile] = tile_normal_tex_paths_[iTile].c_str();
  }

  CreateUniformBuffer(m_pDevice, sizeof(CameraAttribs), "Camera attribs CB", &camera_attribs_);
  CreateUniformBuffer(m_pDevice, sizeof(LightAttribs), "Light attribs CB", &light_attribs_);

  const auto& swap_chain_desc = m_pSwapChain->GetDesc();
  epipolar_light_scattering_.reset(
    new EpipolarLightScattering(m_pDevice, m_pImmediateContext, swap_chain_desc.ColorBufferFormat, swap_chain_desc.DepthBufferFormat, TEX_FORMAT_R11G11B10_FLOAT));
  auto* media_scattering_params = epipolar_light_scattering_->GetMediaAttribsCB();

  earth_hemisphere_.Create(elev_data_source_.get(), terrain_render_params_, m_pDevice, m_pImmediateContext, mtrl_mask_path_.c_str(), tile_diffuse_tex_paths,
                           tile_normal_tex_paths, camera_attribs_, light_attribs_, media_scattering_params);

  CreateShadowMap();
}

void AtmosphereSample::UpdateUI()
{
  ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_FirstUseEver);
  if (ImGui::Begin("设置", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
  {
    ImGui::gizmo3D("光照方向", static_cast<float3&>(light_dir_), ImGui::GetTextLineHeight() * 10);
    ImGui::SliderFloat("相机高度", &camera_pos_.y, 2000, 100000);

    ImGui::SetNextTreeNodeOpen(true, ImGuiCond_FirstUseEver);
    if (ImGui::TreeNode("阴影"))
    {
      {
        constexpr int MinShadowMapSize = 512;
        int           ShadowMapComboId = 0;
        while ((MinShadowMapSize << ShadowMapComboId) != static_cast<int>(shadow_settings_.Resolution)) ++ShadowMapComboId;
        if (ImGui::Combo("图左移数", &ShadowMapComboId,
                         "0\0"
                         "1\0"
                         "2\0"
                         "4\0"
                         "8\0\0"))
        {
          shadow_settings_.Resolution = MinShadowMapSize << ShadowMapComboId;
          CreateShadowMap();
        }
      }

      if (ImGui::SliderInt("级数", &terrain_render_params_.num_shadow_cascades, 1, 8)) CreateShadowMap();

      ImGui::Checkbox("显示分级", &shadow_settings_.visualize_cascades);

      ImGui::TreePop();
    }

    ImGui::Checkbox("启动光散", &enable_light_scattering_);

    if (enable_light_scattering_)
    {
      if (ImGui::BeginTabBar("##tabs", ImGuiTabBarFlags_None))
      {
        if (ImGui::BeginTabItem("基本"))
        {
          ImGui::Checkbox("Enable light shafts", &epipolar_light_scattering_attribs_.bEnableLightShafts);

          static_assert(LIGHT_SCTR_TECHNIQUE_EPIPOLAR_SAMPLING == 0 && LIGHT_SCTR_TECHNIQUE_BRUTE_FORCE == 1, "Unexpected value");
          ImGui::Combo("Light scattering tech", &epipolar_light_scattering_attribs_.iLightSctrTechnique,
                       "Epipolar\0"
                       "Brute force\0\0");

          if (epipolar_light_scattering_attribs_.iLightSctrTechnique == LIGHT_SCTR_TECHNIQUE_EPIPOLAR_SAMPLING)
          {
            {
              static constexpr Uint32 MinSlices    = 128;
              int                     SelectedItem = PlatformMisc::GetLSB(epipolar_light_scattering_attribs_.uiNumEpipolarSlices / MinSlices);
              if (ImGui::Combo("Num Slices", &SelectedItem,
                               "128\0"
                               "256\0"
                               "512\0"
                               "1024\0"
                               "2048\0\0"))
              { epipolar_light_scattering_attribs_.uiNumEpipolarSlices = MinSlices << SelectedItem; }
              ImGui::HelpMarker("Total number of epipolar slices (or lines). For high quality effect, set this value "
                                "to (Screen Width + Screen Height)/2");
            }

            {
              static constexpr Uint32 MinSamples   = 128;
              int                     SelectedItem = PlatformMisc::GetLSB(epipolar_light_scattering_attribs_.uiMaxSamplesInSlice / MinSamples);
              if (ImGui::Combo("Max samples", &SelectedItem,
                               "128\0"
                               "256\0"
                               "512\0"
                               "1024\0"
                               "2048\0\0"))
              { epipolar_light_scattering_attribs_.uiMaxSamplesInSlice = MinSamples << SelectedItem; }
              ImGui::HelpMarker("Maximum number of samples on a single epipolar line. For high quality effect, set "
                                "this value to (Screen Width + Screen Height)/2");
            }

            {
              static constexpr Uint32 MinInitialStep = 4;
              int                     SelectedItem   = PlatformMisc::GetLSB(epipolar_light_scattering_attribs_.uiInitialSampleStepInSlice / MinInitialStep);
              if (ImGui::Combo("Initial Step", &SelectedItem,
                               "4\0"
                               "8\0"
                               "16\0"
                               "32\0"
                               "64\0\0"))
              { epipolar_light_scattering_attribs_.uiInitialSampleStepInSlice = MinInitialStep << SelectedItem; }
              ImGui::HelpMarker("Initial ray marching sample spacing on an epipolar line. Additional samples are added "
                                "at discontinuities.");
            }

            ImGui::SliderFloat("Refinement Threshold", &epipolar_light_scattering_attribs_.fRefinementThreshold, 0.001f, 0.5f);
            ImGui::HelpMarker("Refinement threshold controls detection of discontinuities. Smaller values produce more "
                              "samples and higher quality, but at a higher performance cost.");

            ImGui::Checkbox("Show Sampling", &epipolar_light_scattering_attribs_.bShowSampling);

            if (rg16u_fmt_supported || rg32f_fmt_supported)
            {
              ImGui::Checkbox("Use 1D min/max trees", &epipolar_light_scattering_attribs_.bUse1DMinMaxTree);
              ImGui::HelpMarker("Whether to use 1D min/max binary tree optimization. This improves performance for "
                                "higher shadow map resolution. Test it.");
            }

            ImGui::Checkbox("Optimize Sample Locations", &epipolar_light_scattering_attribs_.bOptimizeSampleLocations);
            ImGui::HelpMarker("Optimize sample locations to avoid oversampling. This should generally be TRUE.");

            ImGui::Checkbox("Correct Scattering At Depth Breaks", &epipolar_light_scattering_attribs_.bCorrectScatteringAtDepthBreaks);
            ImGui::HelpMarker("Whether to correct inscattering at depth discontinuities. Improves quality for additional cost.");

            if (epipolar_light_scattering_attribs_.bCorrectScatteringAtDepthBreaks)
            {
              ImGui::Checkbox("Show Depth Breaks", &epipolar_light_scattering_attribs_.bShowDepthBreaks);
              ImGui::HelpMarker("Whether to display pixels which are classified as depth discontinuities and which will be corrected.");
            }
          }

          ImGui::Checkbox("Lighting Only", &epipolar_light_scattering_attribs_.bShowLightingOnly);

          ImGui::EndTabItem();
        }

        if (ImGui::BeginTabItem("Advanced"))
        {
          if (!epipolar_light_scattering_attribs_.bEnableLightShafts && epipolar_light_scattering_attribs_.iSingleScatteringMode == SINGLE_SCTR_MODE_INTEGRATION)
          {
            ImGui::SliderIntT("Num Integration Steps", &epipolar_light_scattering_attribs_.uiInstrIntegralSteps, 5, 100);
            ImGui::HelpMarker("Number of inscattering integral steps taken when computing unshadowed inscattering");
          }

          if (epipolar_light_scattering_attribs_.iLightSctrTechnique == LIGHT_SCTR_TECHNIQUE_EPIPOLAR_SAMPLING)
          {
            int SelectedItem = PlatformMisc::GetLSB(epipolar_light_scattering_attribs_.uiEpipoleSamplingDensityFactor);
            if (ImGui::Combo("Epipole sampling density", &SelectedItem,
                             "1\0"
                             "2\0"
                             "4\0"
                             "8\0\0"))
            { epipolar_light_scattering_attribs_.uiEpipoleSamplingDensityFactor = 1 << SelectedItem; }
            ImGui::HelpMarker("Sample density scale near the epipole where inscattering changes rapidly. "
                              "Note that sampling near the epipole is very cheap since only a few steps are required "
                              "to perform ray marching.");
          }

          // clang-format off
                    static_assert(SINGLE_SCTR_MODE_NONE         == 0 &&
                                  SINGLE_SCTR_MODE_INTEGRATION  == 1 &&
                                  SINGLE_SCTR_MODE_LUT          == 2, "Unexpected value");
          // clang-format on
          ImGui::Combo("Single scattering mode", &epipolar_light_scattering_attribs_.iSingleScatteringMode,
                       "None\0"
                       "Integration\0"
                       "Look-up table\0\0");

          // clang-format off
                    static_assert(MULTIPLE_SCTR_MODE_NONE        == 0 &&
                                  MULTIPLE_SCTR_MODE_UNOCCLUDED  == 1 &&
                                  MULTIPLE_SCTR_MODE_OCCLUDED    == 2, "Unexpected value");
          // clang-format on
          ImGui::Combo("Higher-order scattering mode", &epipolar_light_scattering_attribs_.iMultipleScatteringMode,
                       "None\0"
                       "Unoccluded\0"
                       "Occluded\0\0");

          // clang-format off
                    static_assert(CASCADE_PROCESSING_MODE_SINGLE_PASS     == 0 &&
                                  CASCADE_PROCESSING_MODE_MULTI_PASS      == 1 &&
                                  CASCADE_PROCESSING_MODE_MULTI_PASS_INST == 2, "Unexpected value");
          // clang-format on
          ImGui::Combo("Cascade processing mode", &epipolar_light_scattering_attribs_.iCascadeProcessingMode,
                       "Single pass\0"
                       "Multi-pass\0"
                       "Multi-pass inst\0\0");

          ImGui::SliderInt("First Cascade to Ray March", &epipolar_light_scattering_attribs_.iFirstCascadeToRayMarch, 0, terrain_render_params_.num_shadow_cascades - 1);
          ImGui::HelpMarker("First cascade to use for ray marching. Usually first few cascades are small, and ray "
                            "marching them is inefficient.");

          if (rg16u_fmt_supported && rg32f_fmt_supported)
          {
            ImGui::Checkbox("32-bit float min/max Shadow Map", &epipolar_light_scattering_attribs_.bIs32BitMinMaxMipMap);
            ImGui::HelpMarker("Whether to use 32-bit float or 16-bit UNORM min-max binary tree.");
          }

          if (epipolar_light_scattering_attribs_.iLightSctrTechnique == LIGHT_SCTR_TECHNIQUE_EPIPOLAR_SAMPLING)
          {
            // clang-format off
                        static_assert(REFINEMENT_CRITERION_DEPTH_DIFF  == 0 &&
                                      REFINEMENT_CRITERION_INSCTR_DIFF == 1, "Unexpected value");
            // clang-format on
            ImGui::Combo("Refinement criterion", &epipolar_light_scattering_attribs_.iRefinementCriterion,
                         "Depth difference\0"
                         "Scattering difference\0\0");
            ImGui::HelpMarker("Epipolar sampling refinement criterion.");

            // clang-format off
                        static_assert(EXTINCTION_EVAL_MODE_PER_PIXEL == 0 &&
                                      EXTINCTION_EVAL_MODE_EPIPOLAR  == 1, "Unexpected value");
            // clang-format on
            ImGui::Combo("Extinction eval mode", &epipolar_light_scattering_attribs_.iExtinctionEvalMode,
                         "Per pixel\0"
                         "Epipolar\0\0");
            ImGui::HelpMarker("Epipolar sampling refinement criterion.");
          }

          if (ImGui::InputFloat("Aerosol Density", &epipolar_light_scattering_attribs_.fAerosolDensityScale, 0.1f, 0.25f, "%.3f", ImGuiInputTextFlags_EnterReturnsTrue))
            epipolar_light_scattering_attribs_.fAerosolDensityScale = clamp(epipolar_light_scattering_attribs_.fAerosolDensityScale, 0.1f, 5.0f);

          if (ImGui::InputFloat("Aerosol Absorption", &epipolar_light_scattering_attribs_.fAerosolAbsorbtionScale, 0.1f, 0.25f, "%.3f",
                                ImGuiInputTextFlags_EnterReturnsTrue))
            epipolar_light_scattering_attribs_.fAerosolAbsorbtionScale = clamp(epipolar_light_scattering_attribs_.fAerosolAbsorbtionScale, 0.0f, 5.0f);

          ImGui::Checkbox("Use custom scattering coeffs", &epipolar_light_scattering_attribs_.bUseCustomSctrCoeffs);
          ImGui::Checkbox("Use Ozone approximation", &epipolar_light_scattering_attribs_.bUseOzoneApproximation);

          if (epipolar_light_scattering_attribs_.bUseCustomSctrCoeffs)
          {
            static constexpr float RLGH_COLOR_SCALE  = 5e-5f;
            static constexpr float MIE_COLOR_SCALE   = 5e-5f;
            static constexpr float OZONE_COLOR_SCALE = 5e-6f;

            {
              float3 RayleighColor = custom_rlgh_beta_ / RLGH_COLOR_SCALE;
              if (ImGui::ColorEdit3("Rayleigh Color", &RayleighColor.r))
              { custom_rlgh_beta_ = max(RayleighColor, float3(1, 1, 1) / 255.f) * RLGH_COLOR_SCALE; }
            }

            {
              float3 MieColor = custom_mie_beta_ / MIE_COLOR_SCALE;
              if (ImGui::ColorEdit3("Mie Color", &MieColor.r)) { custom_mie_beta_ = max(MieColor, float3(1, 1, 1) / 255.f) * MIE_COLOR_SCALE; }
            }

            if (epipolar_light_scattering_attribs_.bUseOzoneApproximation)
            {
              float3 OzoneAbsorption = custom_ozone_absorption_ / OZONE_COLOR_SCALE;
              if (ImGui::ColorEdit3("Ozone Absorbption", &OzoneAbsorption.r))
              { custom_ozone_absorption_ = max(OzoneAbsorption, float3(1, 1, 1) / 255.f) * OZONE_COLOR_SCALE; }
            }

            if (ImGui::Button("Update coefficients"))
            {
              epipolar_light_scattering_attribs_.f4CustomRlghBeta        = custom_rlgh_beta_;
              epipolar_light_scattering_attribs_.f4CustomMieBeta         = custom_mie_beta_;
              epipolar_light_scattering_attribs_.f4CustomOzoneAbsorption = custom_ozone_absorption_;
            }
          }

          ImGui::EndTabItem();
        }

        if (ImGui::BeginTabItem("Tone mapping"))
        {
          {
            std::array<const char*, 7> ToneMappingMode;
            ToneMappingMode[TONE_MAPPING_MODE_EXP]          = "Exp";
            ToneMappingMode[TONE_MAPPING_MODE_REINHARD]     = "Reinhard";
            ToneMappingMode[TONE_MAPPING_MODE_REINHARD_MOD] = "Reinhard Mod";
            ToneMappingMode[TONE_MAPPING_MODE_UNCHARTED2]   = "Uncharted 2";
            ToneMappingMode[TONE_MAPPING_FILMIC_ALU]        = "Filmic ALU";
            ToneMappingMode[TONE_MAPPING_LOGARITHMIC]       = "Logarithmic";
            ToneMappingMode[TONE_MAPPING_ADAPTIVE_LOG]      = "Adaptive log";
            ImGui::Combo("Tone Mapping Mode", &epipolar_light_scattering_attribs_.ToneMapping.iToneMappingMode, ToneMappingMode.data(),
                         static_cast<int>(ToneMappingMode.size()));
          }

          if (epipolar_light_scattering_attribs_.ToneMapping.iToneMappingMode == TONE_MAPPING_MODE_REINHARD_MOD ||
              epipolar_light_scattering_attribs_.ToneMapping.iToneMappingMode == TONE_MAPPING_MODE_UNCHARTED2 ||
              epipolar_light_scattering_attribs_.ToneMapping.iToneMappingMode == TONE_MAPPING_LOGARITHMIC ||
              epipolar_light_scattering_attribs_.ToneMapping.iToneMappingMode == TONE_MAPPING_ADAPTIVE_LOG)
          { ImGui::SliderFloat("White Point", &epipolar_light_scattering_attribs_.ToneMapping.fWhitePoint, 0.01f, 10.0f); }

          if (epipolar_light_scattering_attribs_.ToneMapping.iToneMappingMode == TONE_MAPPING_MODE_EXP ||
              epipolar_light_scattering_attribs_.ToneMapping.iToneMappingMode == TONE_MAPPING_MODE_REINHARD ||
              epipolar_light_scattering_attribs_.ToneMapping.iToneMappingMode == TONE_MAPPING_MODE_REINHARD_MOD ||
              epipolar_light_scattering_attribs_.ToneMapping.iToneMappingMode == TONE_MAPPING_LOGARITHMIC ||
              epipolar_light_scattering_attribs_.ToneMapping.iToneMappingMode == TONE_MAPPING_ADAPTIVE_LOG)
          { ImGui::SliderFloat("Luminance Saturation", &epipolar_light_scattering_attribs_.ToneMapping.fLuminanceSaturation, 0.01f, 2.f); }

          ImGui::SliderFloat("Middle Gray", &epipolar_light_scattering_attribs_.ToneMapping.fMiddleGray, 0.01f, 1.f);
          ImGui::Checkbox("Auto Exposure", &epipolar_light_scattering_attribs_.ToneMapping.bAutoExposure);
          if (epipolar_light_scattering_attribs_.ToneMapping.bAutoExposure) ImGui::Checkbox("Light Adaptation", &epipolar_light_scattering_attribs_.ToneMapping.bLightAdaptation);

          ImGui::EndTabItem();
        }

        ImGui::EndTabBar();
      }
    }
  }
  ImGui::End();
}

AtmosphereSample::~AtmosphereSample() {}

void AtmosphereSample::CreateShadowMap()
{
  ShadowMapManager::InitInfo shadow_map_manager_init_info;
  shadow_map_manager_init_info.Format      = terrain_render_params_.ShadowMapFormat;
  shadow_map_manager_init_info.Resolution  = shadow_settings_.Resolution;
  shadow_map_manager_init_info.NumCascades = terrain_render_params_.num_shadow_cascades;
  shadow_map_manager_init_info.ShadowMode  = SHADOW_MODE_PCF;

  if (!comparison_sampler_)
  {
    SamplerDesc comparsion_sampler_desc;
    comparsion_sampler_desc.ComparisonFunc = COMPARISON_FUNC_LESS;
    // Note: anisotropic filtering requires SampleGrad to fix artifacts at
    // cascade boundaries
    comparsion_sampler_desc.MinFilter = FILTER_TYPE_COMPARISON_LINEAR;
    comparsion_sampler_desc.MagFilter = FILTER_TYPE_COMPARISON_LINEAR;
    comparsion_sampler_desc.MipFilter = FILTER_TYPE_COMPARISON_LINEAR;
    m_pDevice->CreateSampler(comparsion_sampler_desc, &comparison_sampler_);
  }
  shadow_map_manager_init_info.pComparisonSampler = comparison_sampler_;

  shadow_mag_manager_.Initialize(m_pDevice, shadow_map_manager_init_info);
}

void AtmosphereSample::RenderShadowMap(IDeviceContext* pContext, LightAttribs& LightAttribs, const float4x4& mCameraView, const float4x4& mCameraProj)
{
  auto& ShadowAttribs = LightAttribs.ShadowAttribs;

  ShadowMapManager::DistributeCascadeInfo DistrInfo;
  DistrInfo.pCameraView         = &mCameraView;
  DistrInfo.pCameraProj         = &mCameraProj;
  DistrInfo.pLightDir           = &light_dir_;
  DistrInfo.fPartitioningFactor = 0.95f;
  DistrInfo.SnapCascades        = true;
  DistrInfo.EqualizeExtents     = true;
  DistrInfo.StabilizeExtents    = true;
  DistrInfo.AdjustCascadeRange  = [this](int iCascade, float& MinZ, float& MaxZ) {
    if (iCascade < 0)
    {
      // Snap camera z range to the exponential scale
      const float pw = 1.1f;
      MinZ           = std::pow(pw, std::floor(std::log(std::max(MinZ, 1.f)) / std::log(pw)));
      MinZ           = std::max(MinZ, 10.f);
      MaxZ           = std::pow(pw, std::ceil(std::log(std::max(MaxZ, 1.f)) / std::log(pw)));
    }
    else if (iCascade == epipolar_light_scattering_attribs_.iFirstCascadeToRayMarch)
    {
      // Ray marching always starts at the camera position, not at the near plane.
      // So we must make sure that the first cascade used for ray marching covers the camera position
      MinZ = 10.f;
    }
  };
  shadow_mag_manager_.DistributeCascades(DistrInfo, ShadowAttribs);

  // Render cascades
  for (int iCascade = 0; iCascade < terrain_render_params_.num_shadow_cascades; ++iCascade)
  {
    auto* pCascadeDSV = shadow_mag_manager_.GetCascadeDSV(iCascade);

    m_pImmediateContext->SetRenderTargets(0, nullptr, pCascadeDSV, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);
    m_pImmediateContext->ClearDepthStencil(pCascadeDSV, CLEAR_DEPTH_FLAG, 1.f, 0, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);

    const auto CascadeProjMatr = shadow_mag_manager_.GetCascadeTranform(iCascade).Proj;

    auto WorldToLightViewSpaceMatr = ShadowAttribs.mWorldToLightViewT.Transpose();
    auto WorldToLightProjSpaceMatr = WorldToLightViewSpaceMatr * CascadeProjMatr;

    {
      MapHelper<CameraAttribs> CamAttribs(m_pImmediateContext, camera_attribs_, MAP_WRITE, MAP_FLAG_DISCARD);
      CamAttribs->mViewProjT = WorldToLightProjSpaceMatr.Transpose();
    }

    earth_hemisphere_.Render(m_pImmediateContext, terrain_render_params_, WorldToLightProjSpaceMatr, nullptr, nullptr, nullptr, true);
  }
}

// Render a frame
void AtmosphereSample::Render()
{
  float4x4 mViewProj = m_mCameraView * m_mCameraProj;

  LightAttribs LightAttrs;
  LightAttrs.f4Direction   = light_dir_;
  LightAttrs.f4Direction.w = 0;

  float4 f4ExtraterrestrialSunColor = float4(10, 10, 10, 10);
  LightAttrs.f4Intensity            = f4ExtraterrestrialSunColor; // *m_fScatteringScale;
  LightAttrs.f4AmbientLight         = float4(0, 0, 0, 0);

  LightAttrs.ShadowAttribs.iNumCascades = terrain_render_params_.num_shadow_cascades;
  if (shadow_settings_.Resolution >= 2048) LightAttrs.ShadowAttribs.fFixedDepthBias = 0.0025f;
  else if (shadow_settings_.Resolution >= 1024)
    LightAttrs.ShadowAttribs.fFixedDepthBias = 0.0050f;
  else
    LightAttrs.ShadowAttribs.fFixedDepthBias = 0.0075f;

  // m_iFirstCascade must be initialized before calling RenderShadowMap()!
  epipolar_light_scattering_attribs_.iFirstCascadeToRayMarch = std::min(epipolar_light_scattering_attribs_.iFirstCascadeToRayMarch, terrain_render_params_.num_shadow_cascades - 1);

  RenderShadowMap(m_pImmediateContext, LightAttrs, m_mCameraView, m_mCameraProj);

  LightAttrs.ShadowAttribs.bVisualizeCascades = shadow_settings_.visualize_cascades;

  {
    MapHelper<LightAttribs> LightAttribsCBData(m_pImmediateContext, light_attribs_, MAP_WRITE, MAP_FLAG_DISCARD);
    *LightAttribsCBData = LightAttrs;
  }

  // The first time GetAmbientSkyLightSRV() is called, the ambient sky light texture
  // is computed and render target is set. So we need to query the texture before setting
  // render targets
  auto* pAmbientSkyLightSRV = epipolar_light_scattering_->GetAmbientSkyLightSRV(m_pDevice, m_pImmediateContext);

  auto* pRTV = m_pSwapChain->GetCurrentBackBufferRTV();
  auto* pDSV = m_pSwapChain->GetDepthBufferDSV();
  m_pImmediateContext->SetRenderTargets(1, &pRTV, pDSV, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);

  const float ClearColor[] = {0.350f, 0.350f, 0.350f, 1.0f};
  const float Zero[]       = {0.f, 0.f, 0.f, 0.f};
  m_pImmediateContext->ClearRenderTarget(pRTV, enable_light_scattering_ ? Zero : ClearColor, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);

  if (enable_light_scattering_)
  {
    pRTV = m_pOffscreenColorBuffer->GetDefaultView(TEXTURE_VIEW_RENDER_TARGET);
    pDSV = m_pOffscreenDepthBuffer->GetDefaultView(TEXTURE_VIEW_DEPTH_STENCIL);
    m_pImmediateContext->SetRenderTargets(1, &pRTV, pDSV, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);
    m_pImmediateContext->ClearRenderTarget(pRTV, Zero, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);
  }

  m_pImmediateContext->ClearDepthStencil(pDSV, CLEAR_DEPTH_FLAG, 1.f, 0, RESOURCE_STATE_TRANSITION_MODE_TRANSITION);

  CameraAttribs CamAttribs;
  CamAttribs.mViewT        = m_mCameraView.Transpose();
  CamAttribs.mProjT        = m_mCameraProj.Transpose();
  CamAttribs.mViewProjT    = mViewProj.Transpose();
  CamAttribs.mViewProjInvT = mViewProj.Inverse().Transpose();
  float fNearPlane = 0.f, fFarPlane = 0.f;
  m_mCameraProj.GetNearFarClipPlanes(fNearPlane, fFarPlane, is_gl_device_);
  CamAttribs.fNearPlaneZ      = fNearPlane;
  CamAttribs.fFarPlaneZ       = fFarPlane * 0.999999f;
  CamAttribs.f4Position       = camera_pos_;
  CamAttribs.f4ViewportSize.x = static_cast<float>(m_pSwapChain->GetDesc().Width);
  CamAttribs.f4ViewportSize.y = static_cast<float>(m_pSwapChain->GetDesc().Height);
  CamAttribs.f4ViewportSize.z = 1.f / CamAttribs.f4ViewportSize.x;
  CamAttribs.f4ViewportSize.w = 1.f / CamAttribs.f4ViewportSize.y;

  {
    MapHelper<CameraAttribs> CamAttribsCBData(m_pImmediateContext, camera_attribs_, MAP_WRITE, MAP_FLAG_DISCARD);
    *CamAttribsCBData = CamAttribs;
  }

  // Render terrain
  auto* pPrecomputedNetDensitySRV = epipolar_light_scattering_->GetPrecomputedNetDensitySRV();
  terrain_render_params_.DstRTVFormat =
    enable_light_scattering_ ? m_pOffscreenColorBuffer->GetDesc().Format : m_pSwapChain->GetDesc().ColorBufferFormat;
  earth_hemisphere_.Render(m_pImmediateContext, terrain_render_params_, mViewProj, shadow_mag_manager_.GetSRV(), pPrecomputedNetDensitySRV,
                           pAmbientSkyLightSRV, false);

  if (enable_light_scattering_)
  {
    EpipolarLightScattering::FrameAttribs FrameAttribs;

    FrameAttribs.pDevice        = m_pDevice;
    FrameAttribs.pDeviceContext = m_pImmediateContext;
    FrameAttribs.dElapsedTime   = m_fElapsedTime;
    FrameAttribs.pLightAttribs  = &LightAttrs;
    FrameAttribs.pCameraAttribs = &CamAttribs;

    epipolar_light_scattering_attribs_.iNumCascades = terrain_render_params_.num_shadow_cascades;
    epipolar_light_scattering_attribs_.fNumCascades = (float)terrain_render_params_.num_shadow_cascades;

    FrameAttribs.pcbLightAttribs  = light_attribs_;
    FrameAttribs.pcbCameraAttribs = camera_attribs_;

    epipolar_light_scattering_attribs_.fMaxShadowMapStep = static_cast<float>(shadow_settings_.Resolution / 4);

    epipolar_light_scattering_attribs_.f2ShadowMapTexelSize =
      float2(1.f / static_cast<float>(shadow_settings_.Resolution), 1.f / static_cast<float>(shadow_settings_.Resolution));
    epipolar_light_scattering_attribs_.uiMaxSamplesOnTheRay = shadow_settings_.Resolution;

    epipolar_light_scattering_attribs_.uiNumSamplesOnTheRayAtDepthBreak = 32u;

    // During the ray marching, on each step we move by the texel size in either horz
    // or vert direction. So resolution of min/max mipmap should be the same as the
    // resolution of the original shadow map
    epipolar_light_scattering_attribs_.uiMinMaxShadowMapResolution    = shadow_settings_.Resolution;
    epipolar_light_scattering_attribs_.uiInitialSampleStepInSlice     = std::min(epipolar_light_scattering_attribs_.uiInitialSampleStepInSlice, epipolar_light_scattering_attribs_.uiMaxSamplesInSlice);
    epipolar_light_scattering_attribs_.uiEpipoleSamplingDensityFactor = std::min(epipolar_light_scattering_attribs_.uiEpipoleSamplingDensityFactor, epipolar_light_scattering_attribs_.uiInitialSampleStepInSlice);

    FrameAttribs.ptex2DSrcColorBufferSRV = m_pOffscreenColorBuffer->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE);
    FrameAttribs.ptex2DSrcDepthBufferSRV = m_pOffscreenDepthBuffer->GetDefaultView(TEXTURE_VIEW_SHADER_RESOURCE);
    FrameAttribs.ptex2DDstColorBufferRTV = m_pSwapChain->GetCurrentBackBufferRTV();
    FrameAttribs.ptex2DDstDepthBufferDSV = m_pSwapChain->GetDepthBufferDSV();
    FrameAttribs.ptex2DShadowMapSRV      = shadow_mag_manager_.GetSRV();

    // Begin new frame
    epipolar_light_scattering_->PrepareForNewFrame(FrameAttribs, epipolar_light_scattering_attribs_);

    // Render the sun
    epipolar_light_scattering_->RenderSun(pRTV->GetDesc().Format, pDSV->GetDesc().Format, 1);

    // Perform the post processing
    epipolar_light_scattering_->PerformPostProcessing();
  }
}

void GetRaySphereIntersection(float3        f3RayOrigin,
                              const float3& f3RayDirection,
                              const float3& f3SphereCenter,
                              float         fSphereRadius,
                              float2&       f2Intersections)
{
  // http://wiki.cgsociety.org/index.php/Ray_Sphere_Intersection
  f3RayOrigin -= f3SphereCenter;
  float A = dot(f3RayDirection, f3RayDirection);
  float B = 2 * dot(f3RayOrigin, f3RayDirection);
  float C = dot(f3RayOrigin, f3RayOrigin) - fSphereRadius * fSphereRadius;
  float D = B * B - 4 * A * C;
  // If discriminant is negative, there are no real roots hence the ray misses the
  // sphere
  if (D < 0) { f2Intersections = float2(-1, -1); }
  else
  {
    D = sqrt(D);

    f2Intersections = float2(-B - D, -B + D) / (2 * A); // A must be positive here!!
  }
}

void ComputeApproximateNearFarPlaneDist(const float3&   CameraPos,
                                        const float4x4& ViewMatr,
                                        const float4x4& ProjMatr,
                                        const float3&   EarthCenter,
                                        float           fEarthRadius,
                                        float           fMinRadius,
                                        float           fMaxRadius,
                                        float&          fNearPlaneZ,
                                        float&          fFarPlaneZ)
{
  float4x4 ViewProjMatr = ViewMatr * ProjMatr;
  float4x4 ViewProjInv  = ViewProjMatr.Inverse();

  // Compute maximum view distance for the current camera altitude
  float3 f3CameraGlobalPos   = CameraPos - EarthCenter;
  float  fCameraElevationSqr = dot(f3CameraGlobalPos, f3CameraGlobalPos);
  float  fMaxViewDistance    = (float)(sqrt((double)fCameraElevationSqr - (double)fEarthRadius * fEarthRadius) +
                                   sqrt((double)fMaxRadius * fMaxRadius - (double)fEarthRadius * fEarthRadius));
  float  fCameraElev         = sqrt(fCameraElevationSqr);

  fNearPlaneZ = 50.f;
  if (fCameraElev > fMaxRadius)
  {
    // Adjust near clipping plane
    fNearPlaneZ = (fCameraElev - fMaxRadius) / sqrt(1 + 1.f / (ProjMatr._11 * ProjMatr._11) + 1.f / (ProjMatr._22 * ProjMatr._22));
  }

  fNearPlaneZ = std::max(fNearPlaneZ, 50.f);
  fFarPlaneZ  = 1000;

  const int iNumTestDirections = 5;
  for (int i = 0; i < iNumTestDirections; ++i)
  {
    for (int j = 0; j < iNumTestDirections; ++j)
    {
      float3 PosPS, PosWS, DirFromCamera;
      PosPS.x = (float)i / (float)(iNumTestDirections - 1) * 2.f - 1.f;
      PosPS.y = (float)j / (float)(iNumTestDirections - 1) * 2.f - 1.f;
      PosPS.z = 0; // Far plane is at 0 in complimentary depth buffer
      PosWS   = PosPS * ViewProjInv;

      DirFromCamera = PosWS - CameraPos;
      DirFromCamera = normalize(DirFromCamera);

      float2 IsecsWithBottomBoundSphere;
      GetRaySphereIntersection(CameraPos, DirFromCamera, EarthCenter, fMinRadius, IsecsWithBottomBoundSphere);

      float fNearIsecWithBottomSphere = IsecsWithBottomBoundSphere.x > 0 ? IsecsWithBottomBoundSphere.x : IsecsWithBottomBoundSphere.y;
      if (fNearIsecWithBottomSphere > 0)
      {
        // The ray hits the Earth. Use hit point to compute camera space Z
        float3 HitPointWS = CameraPos + DirFromCamera * fNearIsecWithBottomSphere;
        float3 HitPointCamSpace;
        HitPointCamSpace = HitPointWS * ViewMatr;
        fFarPlaneZ       = std::max(fFarPlaneZ, HitPointCamSpace.z);
      }
      else
      {
        // The ray misses the Earth. In that case the whole earth could be seen
        fFarPlaneZ = fMaxViewDistance;
      }
    }
  }
}

void AtmosphereSample::Update(double CurrTime, double ElapsedTime)
{
  const auto& mouseState = m_InputController.GetMouseState();

  float MouseDeltaX = 0;
  float MouseDeltaY = 0;
  if (m_LastMouseState.PosX >= 0 && m_LastMouseState.PosY >= 0 && m_LastMouseState.ButtonFlags != MouseState::BUTTON_FLAG_NONE)
  {
    MouseDeltaX = mouseState.PosX - m_LastMouseState.PosX;
    MouseDeltaY = mouseState.PosY - m_LastMouseState.PosY;
  }
  m_LastMouseState = mouseState;

  if ((m_LastMouseState.ButtonFlags & MouseState::BUTTON_FLAG_LEFT) != 0)
  {
    constexpr float CameraRotationSpeed = 0.005f;
    m_fCameraYaw += MouseDeltaX * CameraRotationSpeed;
    m_fCameraPitch += MouseDeltaY * CameraRotationSpeed;
  }
  m_CameraRotation =
    Quaternion::RotationFromAxisAngle(float3{1, 0, 0}, -m_fCameraPitch) * Quaternion::RotationFromAxisAngle(float3{0, 1, 0}, -m_fCameraYaw);
  camera_pos_.y += mouseState.WheelDelta * 500.f;
  camera_pos_.y = std::max(camera_pos_.y, 2000.f);
  camera_pos_.y = std::min(camera_pos_.y, 100000.f);

  auto CameraRotationMatrix = m_CameraRotation.ToMatrix();

  if ((m_LastMouseState.ButtonFlags & MouseState::BUTTON_FLAG_RIGHT) != 0)
  {
    constexpr float LightRotationSpeed = 0.001f;

    float  fYawDelta   = MouseDeltaX * LightRotationSpeed;
    float  fPitchDelta = MouseDeltaY * LightRotationSpeed;
    float3 WorldUp{CameraRotationMatrix._12, CameraRotationMatrix._22, CameraRotationMatrix._32};
    float3 WorldRight{CameraRotationMatrix._11, CameraRotationMatrix._21, CameraRotationMatrix._31};
    light_dir_ = float4(light_dir_, 0) * float4x4::RotationArbitrary(WorldUp, fYawDelta) * float4x4::RotationArbitrary(WorldRight, fPitchDelta);
  }

  SampleBase::Update(CurrTime, ElapsedTime);
  UpdateUI();

  m_fElapsedTime = static_cast<float>(ElapsedTime);

  const auto& SCDesc = m_pSwapChain->GetDesc();
  // Set world/view/proj matrices and global shader constants
  float aspectRatio = (float)SCDesc.Width / SCDesc.Height;

  m_mCameraView = float4x4::Translation(-camera_pos_) * CameraRotationMatrix;

  // This projection matrix is only used to set up directions in view frustum
  // Actual near and far planes are ignored
  float    FOV      = PI_F / 4.f;
  float4x4 mTmpProj = float4x4::Projection(FOV, aspectRatio, 50.f, 500000.f, is_gl_device_);

  float  fEarthRadius = AirScatteringAttribs().fEarthRadius;
  float3 EarthCenter(0, -fEarthRadius, 0);
  float  fNearPlaneZ, fFarPlaneZ;
  ComputeApproximateNearFarPlaneDist(camera_pos_, m_mCameraView, mTmpProj, EarthCenter, fEarthRadius, fEarthRadius + min_elevation_,
                                     fEarthRadius + max_elevation_, fNearPlaneZ, fFarPlaneZ);
  fNearPlaneZ = std::max(fNearPlaneZ, 50.f);
  fFarPlaneZ  = std::max(fFarPlaneZ, fNearPlaneZ + 100.f);
  fFarPlaneZ  = std::max(fFarPlaneZ, 1000.f);

  m_mCameraProj = float4x4::Projection(FOV, aspectRatio, fNearPlaneZ, fFarPlaneZ, is_gl_device_);

#if 0
    if( m_bAnimateSun )
    {
        auto &LightOrientationMatrix = *m_pDirLightOrienationCamera->GetParentMatrix();
        float3 RotationAxis( 0.5f, 0.3f, 0.0f );
        float3 LightDir = m_pDirLightOrienationCamera->GetLook() * -1;
        float fRotationScaler = ( LightDir.y > +0.2f ) ? 50.f : 1.f;
        float4x4 RotationMatrix = float4x4RotationAxis(RotationAxis, 0.02f * (float)deltaSeconds * fRotationScaler);
        LightOrientationMatrix = LightOrientationMatrix * RotationMatrix;
        m_pDirLightOrienationCamera->SetParentMatrix(LightOrientationMatrix);
    }

    float dt = (float)ElapsedTime;
    if (m_Animate && dt > 0 && dt < 0.2f)
    {
        float3 axis;
        float angle = 0;
        AxisAngleFromRotation(axis, angle, m_Rotation);
        if (length(axis) < 1.0e-6f) 
            axis[1] = 1;
        angle += m_AnimationSpeed * dt;
        if (angle >= 2.0f*FLOAT_PI)
            angle -= 2.0f*FLOAT_PI;
        else if (angle <= 0)
            angle += 2.0f*FLOAT_PI;
        m_Rotation = RotationFromAxisAngle(axis, angle);
    }
#endif
}

void AtmosphereSample::WindowResize(Uint32 Width, Uint32 Height)
{
  epipolar_light_scattering_->OnWindowResize(m_pDevice, Width, Height);
  // Flush is required because Intel driver does not release resources until
  // command buffer is flushed. When window is resized, WindowResize() is called for
  // every intermediate window size, and light scattering object creates resources
  // for the new size. This resources are then released by the light scattering object, but
  // not by Intel driver, which results in memory exhaustion.
  m_pImmediateContext->Flush();

  m_pOffscreenColorBuffer.Release();
  m_pOffscreenDepthBuffer.Release();

  TextureDesc ColorBuffDesc;
  ColorBuffDesc.Name      = "Offscreen color buffer";
  ColorBuffDesc.Type      = RESOURCE_DIM_TEX_2D;
  ColorBuffDesc.Width     = Width;
  ColorBuffDesc.Height    = Height;
  ColorBuffDesc.MipLevels = 1;
  ColorBuffDesc.Format    = TEX_FORMAT_R11G11B10_FLOAT;
  ColorBuffDesc.BindFlags = BIND_SHADER_RESOURCE | BIND_RENDER_TARGET;
  m_pDevice->CreateTexture(ColorBuffDesc, nullptr, &m_pOffscreenColorBuffer);

  TextureDesc DepthBuffDesc = ColorBuffDesc;
  DepthBuffDesc.Name        = "Offscreen depth buffer";
  DepthBuffDesc.Format      = TEX_FORMAT_D32_FLOAT;
  DepthBuffDesc.BindFlags   = BIND_SHADER_RESOURCE | BIND_DEPTH_STENCIL;
  m_pDevice->CreateTexture(DepthBuffDesc, nullptr, &m_pOffscreenDepthBuffer);
}

} // namespace Diligent
