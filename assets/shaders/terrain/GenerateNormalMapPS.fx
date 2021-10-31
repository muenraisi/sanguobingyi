
#include "HostSharedTerrainStructs.fxh"
#include "TerrainShadersCommon.fxh"

Texture2D< uint > g_tex2DElevationMap;

cbuffer cbNMGenerationAttribs
{
	NMGenerationAttribs g_normal_generation_attribs;
};


float3 ComputeNormal(int2 i2ElevMapIJ, float f_sample_spacing_interval, int i_mip_level)
{
	int i_mip_width, i_mip_height;
	// This version of GetDimensions() does not work on D3D12. Looks like a bug in shader compiler
	// g_tex2DElevationMap.GetDimensions( i_mip_level, i_mip_width, i_mip_height, Levels );
	g_tex2DElevationMap.GetDimensions(i_mip_width, i_mip_height);
	i_mip_width = i_mip_width >> i_mip_level;
	i_mip_height = i_mip_height >> i_mip_level;

	int i = i2ElevMapIJ.x;
	int i_plus = min(i + 1, i_mip_width - 1);
	int i_minus = max(i - 1, 0);

	int j = i2ElevMapIJ.y;
	int j_plus = min(j + 1, i_mip_height - 1);
	int j_minus = max(j - 1, 0);

# define GET_ELEV(x,y) float( g_tex2DElevationMap.Load(int3(x,y, i_mip_level)) )

#if 1
	float f_height_00 = GET_ELEV(i_minus, j_minus);
	float f_height_10 = GET_ELEV(i, j_minus);
	float f_height_20 = GET_ELEV(i_plus, j_minus);

	float f_height_01 = GET_ELEV(i_minus, j);
	//float f_height_11 = GET_ELEV(  i0, j0 );
	float f_height_21 = GET_ELEV(i_plus, j);

	float f_height_02 = GET_ELEV(i_minus, j_plus);
	float f_height_12 = GET_ELEV(i, j_plus);
	float f_height_22 = GET_ELEV(i_plus, j_plus);

	float3 f3_grad;
	f3_grad.x = (f_height_00 + f_height_01 + f_height_02) - (f_height_20 + f_height_21 + f_height_22);
	f3_grad.y = (f_height_00 + f_height_10 + f_height_20) - (f_height_02 + f_height_12 + f_height_22);
	f3_grad.z = f_sample_spacing_interval * 6.0;
	//f3_grad.x = (3*f_height_00+10*f_height_01+3*f_height_02) - (3*f_height_20+10*f_height_21+3*f_height_22);
	//f3_grad.y = (3*f_height_00+10*f_height_10+3*f_height_20) - (3*f_height_02+10*f_height_12+3*f_height_22);
	//f3_grad.z = f_sample_spacing_interval * 32.f;
#else
	float f_height_1 = GET_ELEV(i_plus, j);
	float f_height_2 = GET_ELEV(i_minus, j);
	float f_height_3 = GET_ELEV(i, j_plus);
	float f_height_4 = GET_ELEV(i, j_minus);

	float3 f3_grad;
	f3_grad.x = f_height_2 - f_height_1;
	f3_grad.y = f_height_4 - f_height_3;
	f3_grad.z = f_sample_spacing_interval * 2.0;
#endif
	f3_grad.xy *= g_normal_generation_attribs.height_scale;
	float3 f3_normal = normalize(f3_grad);

	return f3_normal;
}


void GenerateNormalMapPS(in float4 f4Pos : SV_Position,	out float2 f2_out_normal_xy : SV_Target)
{
	float3 f3_normal = ComputeNormal(int2(f4Pos.xy), g_normal_generation_attribs.m_fSampleSpacingInterval * exp2(float(g_normal_generation_attribs.m_iMIPLevel)), g_normal_generation_attribs.m_iMIPLevel);
	// Only xy components are stored. z component is calculated in the shader
	f2_out_normal_xy = f3_normal.xy * float2(0.5, 0.5) + float2(0.5, 0.5);
}
