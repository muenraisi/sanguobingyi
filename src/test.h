#include <algorithm>


enum class QUAD_TRI_TYPE
{
    UNDEFINED = 0,
    // 01      11
    //  *------*
    //  |   .' |
    //  | .'   |
    //  * -----*
    // 00      10
    DIAG_00_11,
    // 01      11
    //  *------*
    //  | '.   |
    //  |   '. |
    //  * -----*
    // 00      10
    DIAG_01_10
};