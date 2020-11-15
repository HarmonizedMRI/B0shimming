function dur = roundToRaster(dur, raster)

dur = dur + raster - mod(dur, raster);

