

# Our "standard" variable names, and how they map to hycom names
atmosphere_variable_names = {
      "10u"     : "wndewd",
      "10v"     : "wndnwd",
      "2t"      : "airtmp",
      "vapmix"  : "vapmix",
      "msl"     : "mslprs",
      "tp"      : "precip",
      "ssrd"    : "shwflx"
      }


# Units that are required by hycom (udunits)
atmosphere_variable_units = {
      "10u"     : "m s**-1",
      "10v"     : "m s**-1",
      "2t"      : "C",
      "vapmix"  : "kg kj**-1",
      "msl"     : "Pa",
      "tp"      : "m s**-1",
      "ssrd"    : "W m**-2 s*-1"
      }
