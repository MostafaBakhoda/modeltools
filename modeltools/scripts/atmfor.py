#!/usr/bin/env python
import modeltools.hycom
import argparse
import datetime

if __name__ == "__main__" :
   class DateTimeParseAction(argparse.Action) :
       def __call__(self, parser, args, values, option_string=None):
          tmp = datetime.datetime.strptime( values, "%Y-%m-%dT%H:%M:%S")
          setattr(args, self.dest, tmp)
   parser = argparse.ArgumentParser(description='Prepare HYCOM forcing files')
   parser.add_argument('--plot_diag', action="store_true")
   parser.add_argument('start_time', action=DateTimeParseAction, help='Start time in UTC zone. Format = YYYY-mm-ddTHH:MM:SS')
   parser.add_argument('end_time',   action=DateTimeParseAction, help='Stop  time in UTC zone. Format = YYYY-mm-ddTHH:MM:SS')
   parser.add_argument('xml_file',   type=str, help='xml file')
   parser.add_argument('xml_id',   type=str, help='xml id')

   args = parser.parse_args()


   # Set up AtmosphericForcing object, which keeps track of data to be read
   af=modeltools.forcing.atmosphere.AtmosphericForcing(args.xml_file,args.xml_id)



   modeltools.hycom.atmfor(args.start_time,args.end_time,af,plot_diag=args.plot_diag)#,gridfile="regional.grid",blkdat_file="blkdat.input")
