#pragma once

#include "th_raster.h"

namespace WHU {
    enum class CoverType {
        Empty = 0,
    };
    
    using  CoverRaster = GenericRaster<unsigned char>;
}