
#include "TerrainShadersCommon.fxh"
 
void GenerateScreenSizeQuadVS(in uint ui_vertex_id : SV_VertexID,
                              out float4 f4_pos : SV_Position)
{
    float4 f4_min_max_uv = float4(-1.0, -1.0, 1.0, 1.0);
    
    float2 f2_verts[4];
    f2_verts[0] = f4_min_max_uv.xy;
    f2_verts[1] = f4_min_max_uv.xw;
    f2_verts[2] = f4_min_max_uv.zy;
    f2_verts[3] = f4_min_max_uv.zw;

    f4_pos = float4(f2_verts[ui_vertex_id], 1.0, 1.0);
}
