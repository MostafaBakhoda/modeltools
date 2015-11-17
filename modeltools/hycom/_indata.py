import cfunits

# Our "standard" variable names, and how they map to hycom names
atmosphere_variable_names = {
      "10u"     : "wndewd",
      "10v"     : "wndnwd",
      "2t"      : "airtmp",
      "vapmix"  : "vapmix",
      "msl"     : "mslprs",
      "tp"      : "precip",
      "ssrd"    : "shwflx",
      "tsrd"    : "radflx",
      "taux"    : "taux",
      "tauy"    : "tauy"
      }


# Units that are required by hycom (udunits)
atmosphere_variable_units = {
      "10u"     : "m s**-1",
      "10v"     : "m s**-1",
      "2t"      : "degree_Celsius",
      "vapmix"  : "kg kg**-1",
      "msl"     : "Pa",
      "tp"      : "m s**-1",
      "ssrd"    : "W m**-2",
      "tsrd"    : "W m**-2",
      "taux"    : "N m**-2",
      "tauy"    : "N m**-2"
      }


atmosphere_variable_vectors = {
      "taux" : ("taux","tauy"),
      "10u"  : ("10u" ,"10v" )
      }



atmosphere_variable_cfunit = dict(
      [(elem[0],cfunits.Units(elem[1])) for elem in atmosphere_variable_units.items() ]
      )
